import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('TkAgg')
import scienceplots

class Distribution:
    """
    Distribution of MB Energies.
    """
    def __init__(self,path:str) -> None:
        """
        Initalizes a Distribution object by reading a file written by MBAnalysis2.py:
        LifetimesIO.write_counts().

        Parameters
        ----------
        path : str
            Path to the Distribution File.
        """
        data = np.loadtxt(path)
        self.e_mb = data[:,0]
        self.counts = data[:,1]
        bin_width = self.e_mb[1] - self.e_mb[0]
        self.norm_counts = self.counts / (np.sum(self.counts) * bin_width)

    @staticmethod
    def _gauss(x,mu,sigma):
        return (1.0 / (sigma*np.sqrt(2*np.pi))) * np.exp(-0.5*((x-mu) / sigma)**2)
    
    @staticmethod
    def _log_gauss(x,mu,sigma):
        return -np.log(sigma) - 0.5*np.log(2*np.pi) - 0.5*((x-mu) / sigma)**2

    def fit_gauss(self):
        """
        Fits a simple Gauss-Function to the Distribution Data. 
        """
        mu0 = np.sum(self.e_mb * self.counts) / np.sum(self.counts)
        sigma0 = np.sqrt(np.sum(self.counts * (self.e_mb - mu0)**2) / np.sum(self.counts))
        popt, pcov = curve_fit(self._gauss,self.e_mb,self.norm_counts,p0=[mu0,sigma0])
        perr = np.sqrt(np.diag(pcov))

        self.mu = popt[0]
        self.sigma = popt[1]

        return self.mu, self.sigma
    
    def fit_gauss_ln_right(self):
        """
        Fits a logarithmic Gauss-Function to the right flank of the Distribution Data. 
        """
        peak_idx = np.argmax(self.norm_counts)
        e_right = self.e_mb[peak_idx:]
        norm_right = self.norm_counts[peak_idx:]

        mask = norm_right > 0
        e_right = e_right[mask]
        log_norm_right = np.log(norm_right[mask])

        mu0 = self.e_mb[peak_idx]
        sigma0 = (self.e_mb[-1]-mu0) / 3.0

        popt, pcov = curve_fit(self._log_gauss,e_right,log_norm_right,p0=[mu0,sigma0])
        perr = np.sqrt(np.diag(pcov))

        self.mu = popt[0]
        self.sigma = popt[1]

        return self.mu, self.sigma
    
    def fit_gauss_ln(self):
        """
        Fits a logarithmic Gauss-Function to the Distribution Data. 
        """
        mask = self.norm_counts > 0
        e_mb = self.e_mb[mask]
        log_norm_right = np.log(self.norm_counts[mask])

        mu0 = np.sum(self.e_mb * self.counts) / np.sum(self.counts)
        sigma0 = np.sqrt(np.sum(self.counts * (self.e_mb - mu0)**2) / np.sum(self.counts))

        popt, pcov = curve_fit(self._log_gauss,e_mb,log_norm_right,p0=[mu0,sigma0])
        perr = np.sqrt(np.diag(pcov))

        self.mu = popt[0]
        self.sigma = popt[1]

        return self.mu, self.sigma
    
    def fit_gauss_right(self):
        """
        Fits a simple Gauss-Function to the right flank of the Distribution Data. 
        """
        peak_idx = np.argmax(self.norm_counts)
        e_right = self.e_mb[peak_idx:]
        norm_right = self.norm_counts[peak_idx:]

        mu0 = self.e_mb[peak_idx]
        sigma0 = (self.e_mb[-1]-mu0) / 3.0

        popt, pcov = curve_fit(self._gauss,e_right,norm_right,p0=[mu0,sigma0])
        perr = np.sqrt(np.diag(pcov))

        self.mu = popt[0]
        self.sigma = popt[1]

        return self.mu, self.sigma
    
    def plot_gauss(self,path:str="gauss") -> None:
        """
        Plots the underlying Distribution and the used Gaussian Fit.

        Parameters
        ----------
        path : str
            Path to the Output File.
        """
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'sans-serif'
        bin_width = self.e_mb[1] - self.e_mb[0]
        x_fine = np.linspace(self.e_mb[0]-0.01, self.e_mb[-1]+0.01, 500)
        fig, ax = plt.subplots(figsize=(8,5))                          # ← Fix
        ax.scatter(self.e_mb, self.norm_counts, alpha=0.5)
        ax.plot(x_fine, self._gauss(x_fine, self.mu, self.sigma), 'r-')
        ax.set_xlabel(r"$e_{MB}$")
        ax.set_ylabel(r"$G(e_{MB})$")
        plt.tight_layout()
        plt.show()

    def set_temp(self,temp_min:float=0.08,temp_max:float=0.15,temp_step:float=0.005) -> None:
        """
        Sets a temperature range to analyze.

        Parameters
        ----------
        temp_min : float
            Minimal Temperature.
        temp_max : float
            Maximal Temperature.
        temp_step : float
            Temperature Step.
        """
        self.temp = np.arange(temp_min,temp_max,temp_step)
        return
    
    def calc_avg(self,temp:float=0.1) -> float:
        """
        Calculates the average energy <E> for a given temperature.

        Parameters
        ----------
        temp : float
            Temperature to calculate.
        """
        inv_temp = 1/temp
        mask = self.norm_counts > 0
        e = self.e_mb[mask]
        g = self.norm_counts[mask]
        ener_avg = 0.0
        part_sum = 0.0
        for i,eis in enumerate(e):
            exponent = np.log(g[i])-eis*inv_temp
            ener_avg += eis*np.exp(exponent)
            part_sum += np.exp(exponent)
        return (ener_avg/part_sum)
    
    def plot_avg(self,path:str="avg") -> None:
        """
        Plots the average MB energy calculated from the used Gaussian Fit.

        Parameters
        ----------
        path : str
            Path to the Output File.
        """
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'sans-serif'

        fig, ax = plt.subplots(figsize=(8,5)) 
        x = [(1/temp) for i,temp in enumerate(self.temp) if np.isfinite(self.calc_avg(temp=temp))]
        y = [self.calc_avg(temp=temp) for i,temp in enumerate(self.temp)]

        ax.scatter(x,y)
        ax.set_ylabel(r"$\langle E\rangle$")
        ax.set_xlabel(r"$1/T$")
        plt.tight_layout()
        plt.show()


dst = Distribution(path="04_Comp/sgk_N66-IS.dat")
dst.fit_gauss()
dst.plot_gauss()
dst.set_temp()
dst.plot_avg()