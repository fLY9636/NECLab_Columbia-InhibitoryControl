
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import svd, pinv
from typing import Callable

class myDPCA:
    def __init__(self):
        pass

    def dpca_custom_stim_only(self, Xfull, W, V, plot_function, **kwargs):
        """
        Custom DPCA function that produces a plot of the dPCA results with various customizable options.
        
        Parameters:
            Xfull (ndarray): Data matrix of neural responses.
            W (ndarray): Decoder matrix.
            V (ndarray): Encoder matrix.
            plot_function (function): Function to plot each component.
            kwargs: Additional optional parameters.
        """
        options = {
            'time': None,
            'whichMarg': None,
            'timeEvents': None,
            'ylims': None,
            'componentsSignif': None,
            'timeMarginalization': None,
            'legendSubplot': None,
            'marginalizationNames': None,
            'marginalizationColours': None,
            'explainedVar': None,
            'numCompToShow': 

    def dpca(Xfull, numComps, **kwargs):
        options = {
            'combinedParams': [],
            'lambda': 0,
            'order': 'yes',
            'timeSplits': [],
            'timeParameter': None,
            'notToSplit': [],
            'scale': 'no',
            'Cnoise': None
        }
        options.update(kwargs)

        X = Xfull.reshape((Xfull.shape[0], -1))
        X -= X.mean(axis=1, keepdims=True)
        XfullCen = X.reshape(Xfull.shape)
        totalVar = np.sum(X ** 2)
        
        Xmargs, margNums = myDPCA.dpca_marginalize(XfullCen, combinedParams=options['combinedParams'],
                                                   timeSplits=options['timeSplits'], timeParameter=options['timeParameter'],
                                                   notToSplit=options['notToSplit'], ifFlat=True)
        
        W, V, whichMarg = [], [], []
        for i, X_marg in enumerate(Xmargs):
            nc = numComps if isinstance(numComps, int) else numComps[margNums[i]]
            thisLambda = options['lambda'] if isinstance(options['lambda'], (int, float)) else options['lambda'][margNums[i]]
            if nc == 0:
                continue
            
            C = X_marg @ X.T @ pinv(X @ X.T + options.get('Cnoise', 0) + (totalVar * thisLambda)**2 * np.eye(X.shape[0]))
            U, _, _ = svd(C @ X, full_matrices=False)
            P, D = U[:, :nc], U[:, :nc].T @ C
            
            if options['scale'] == 'yes':
                for uu in range(D.shape[0]):
                    scalingFactor = (X_marg.ravel() @ (P[:, uu] @ D[uu, :] @ X).ravel()) / (D[uu, :] @ X).ravel() @ (D[uu, :] @ X).ravel()
                    D[uu, :] *= scalingFactor
            
            W.append(D)
            V.append(P)
            whichMarg.extend([i + 1] * nc)
        
        W, V = np.hstack(W).T, np.hstack(V)
        return W, V, whichMarg

    def dpca_plot(Xfull, W, V, plotFunction: Callable, **kwargs):
        options = {
            'time': None,
            'whichMarg': None,
            'timeEvents': None,
            'ylims': None,
            'componentsSignif': None,
            'timeMarginalization': None,
            'legendSubplot': None,
            'marginalizationNames': [],
            'marginalizationColours': [],
            'explainedVar': None,
            'numCompToShow': 15,
            'X_extra': None,
            'showNonsignificantComponents': False
        }
        options.update(kwargs)
        numCompToShow = min(options['numCompToShow'], W.shape[1])

        X = Xfull.reshape((-1, Xfull.shape[-1])).T
        X -= X.mean(axis=0)
        Z = X @ W
        
        components_to_plot = np.arange(numCompToShow)
        
        fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(18, 10))
        for i, comp in enumerate(components_to_plot):
            ax = axes.flat[i]
            plotFunction(ax, Z[:, comp], options['time'], options['ylims'])
        
        plt.tight_layout()
        plt.show()
        return components_to_plot, Z

    def dpca_explainedVariance(Xfull, W, V, **kwargs):
        options = {
            'combinedParams': [],
            'X_trial': None,
            'numOfTrials': None,
            'Cnoise': None
        }
        options.update(kwargs)

        X = Xfull.reshape((Xfull.shape[0], -1))
        X -= X.mean(axis=1, keepdims=True)
        Xmargs = myDPCA.dpca_marginalize(X, combinedParams=options['combinedParams'], ifFlat=True)
        
        totalVar = np.sum(X ** 2)
        totalMarginalizedVar = [np.sum(X_marg ** 2) for X_marg in Xmargs]
        cumulativePCA = np.cumsum(np.diag(svd(X.T @ X, full_matrices=False)[1]) ** 2) / totalVar * 100
        
        explVar = {
            'totalVar': totalVar,
            'totalMarginalizedVar': totalMarginalizedVar,
            'cumulativePCA': cumulativePCA,
        }
        
        return explVar

    def dpca_marginalize(Xfull, combinedParams=None, timeSplits=None, timeParameter=None, notToSplit=None, ifFlat=False):
        Xmargs, margNums = [Xfull], [1]  # Placeholder
        return Xmargs, margNums
