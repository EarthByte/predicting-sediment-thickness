""" Geo Analysis Preprocessing."""
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.model_selection._split import _BaseKFold

class LongFold(_BaseKFold):
    """
        Design special k-fold: LongFold
        Split X along longitude
    """
    # Initialization
    def __init__(self, n_splits=3, shuffle=False, random_state=None):
        super(LongFold, self).__init__(n_splits, shuffle, random_state)

    def _iter_test_indices(self, X, y=None, groups=None):
        """
        Return test and train indice

        Parameters
        ----------
        X : ndarray (n, d)
            First column must be longitude
        """
        n_clusters = self.n_splits
        n_samples = X.shape[0]

        # original indices
        indices = np.arange(n_samples)

        # extract longitude from X
        lon = X[:, 0]

        # get quantile of longitude given n_splits
        unique_lon = np.unique(lon)
        unique_lon = np.array_split(unique_lon, n_clusters)
        quantiles = []
        for u_lon in unique_lon:
            quantiles.append(u_lon[-1])

        # concatenate lon with indices
        lon_indices = np.column_stack([lon, indices])

        # sort lon_indices based on lon
        sort_lon_indices = lon_indices[lon_indices[:, 0].argsort()]
        sort_lon_indices = sort_lon_indices.astype(int)

        # yield test index
        last_index = np.array([])
        for q in quantiles:
            index = np.where(sort_lon_indices[:, 0] <= q)
            yield np.setdiff1d(sort_lon_indices[index, 1], last_index)
            last_index = sort_lon_indices[index, 1]


def regression(data, regressor, n_splits=3,
               lon_ind=0, lat_ind=1, y_ind=-1,
               fit_with_lon_lat=False, logy=True):
    """
    Conduct regression on data given regressor based on LongFold cross validation.

    Parameters
    ----------
    data : ndarray (n, d)
           Include all features and the target variable
           It is recommended to have longitude in the first column

    regressor : a model

    n_splits : the parameter for LongFold; Optional

    lon_ind : int; column index of longitude, recommended 0; Optional

    lat_ind : int; column index of latitude, recommended 1; Optional

    y_ind : int; column index of target variable, recommend the last column;
            Optional

    fit_with lon_lat : boolean;
                       --if true, then add lon and lat as two features
                       --if false, not add lon and lat as two features

    logy : boolean: if true then log-transform target variable y
    """
    # feature maxtrix - X
    X = data[:, :y_ind]

    # if logy is true, then take log transformation on target variable
    if logy:
        y = np.log(data[:, y_ind])
    # otherwise, just extract target variable from data
    else:
        y = data[:, y_ind]

    # obtain training and testing indices
    kf = LongFold(n_splits=n_splits)
    train_indices = []
    test_indices = []
    for train_index, test_index in kf.split(X):
        train_indices.append(train_index)
        test_indices.append(test_index)

    # if fit_with_lon_lat is false, then delete longitude and latitude columns
    # from feature matrix X
    if ~fit_with_lon_lat:
        X = scipy.delete(X, [lon_ind, lat_ind], 1)

    # initialize predicted y
    y_pred_whole = np.empty([y.shape[0], ])

    # get training and testing indices then do machine learning
    for i in range(len(train_indices)):
        train_index = train_indices[i].tolist()
        test_index = test_indices[i].tolist()

        X_train = X[train_index, :]
        y_train = y[train_index]
        X_test = X[test_index, :]

        # machine learning
        regressor.fit(X_train, y_train)

        # obtain prediction for this test indice
        y_pred = regressor.predict(X_test)

        # fill in y_pred_whole
        y_pred_whole[test_index] = y_pred

    # return predicted y as a whole and trained regressor
    return y_pred_whole, regressor


def predict(X, regressor, lon_ind=0, lat_ind=1, n_splits=3):

    # obtain training and testing indices
    kf = LongFold(n_splits=n_splits)
    train_indices = []
    test_indices = []
    for train_index, test_index in kf.split(X):
        train_indices.append(train_index)
        test_indices.append(test_index)

    X = scipy.delete(X, [lon_ind, lat_ind], 1)

    # machine learning with feature and obtain predicted log thickness
    y_pred_whole = np.empty([X.shape[0], ])

    # get training and testing indices then do machine learning
    for i in range(len(test_indices)):
        test_index = test_indices[i].tolist()

        X_test = X[test_index, :]

        # predict
        y_pred = regressor.predict(X_test)

        y_pred_whole[test_index] = y_pred

    return y_pred_whole


def feature_scatter(y_actual, y_pred, title, xlim1=2, xlim2=10, ylim1=2, ylim2=10):
    """
    Plot the y_pred vs y_actual scatterplot and print the R squared.
    """
    f, ax = plt.subplots()
    ax.scatter(y_actual, y_pred, alpha=0.01, color='r')
    ax.set_xlabel('actual y')
    ax.set_ylabel('predicted y')
    lim1 = xlim1
    lim2 = xlim2
    if ylim1 < lim1:
        lim1 = ylim1
    if ylim2 > lim2:
        lim2 = ylim2
    ax.set_xlim([lim1, lim2])
    ax.set_ylim([lim1, lim2])
    ax.set_title(title + ' with R sqaured = ' + str(r2_score(y_actual, y_pred)))
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", color='r')
    plt.show()


def draw_global(c, title, longitude, latitude, vmin=2.4, vmax=9,
                cmap='RdBu_r'):
    """Draw global map of actual thickness."""
    sc = plt.scatter(longitude, latitude, c=c, cmap=cmap, lw=0,
                     vmin=vmin, vmax=vmax)
    plt.ylim([-70, 78])
    plt.xlim([-180, 180])
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.title(title)
    plt.colorbar(sc)
    plt.show()


def two_feature_analysis(data, regressor, feature1_ind, feature2_ind,
                         feature1_name, feature2_name, y_name,
                         lon_ind=0, lat_ind=1, y_ind=-1, query_size=10,
                         ):
    # set up variables
    X = data[:, :y_ind]
    log_thick = np.log(data[:, y_ind])

    # get feature column
    feature = X[:, [feature1_ind, feature2_ind]]

    # learn y just based on feature
    model = regressor.fit(feature, log_thick)
    print('Polynomial feature names: {0}'.format(model.get_params()['poly'].get_feature_names()))
    print('Polynomial coefficients: {0}'.format(model.get_params()['linear'].coef_))

    # generate query sets
    fea1_min = np.min(feature[:, 0])
    fea1_max = np.max(feature[:, 0])
    fea1_query = np.linspace(fea1_min, fea1_max, query_size)[:, np.newaxis]
    fea2_min = np.min(feature[:, 1])
    fea2_max = np.max(feature[:, 1])
    fea2_query = np.linspace(fea2_min, fea2_max, query_size)[:, np.newaxis]

    # draw regression surface
    xeval, yeval = np.meshgrid(fea1_query, fea2_query)
    query = np.vstack((xeval.ravel(), yeval.ravel())).T
    fea_pred = model.predict(query)
    fea_pred_arr = fea_pred.reshape(query_size, query_size)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xeval, yeval, fea_pred_arr, rstride=1, cstride=1,
                           cmap='coolwarm', lw=0.1, vmin=-4, vmax=4)
    ax.set_xlabel(feature1_name)
    ax.set_ylabel(feature2_name)
    ax.set_zlabel(y_name)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title('Regression Surface')
    fig.savefig('3d_plot.png',bbox_inches='tight',dpi=600)
    plt.show()

