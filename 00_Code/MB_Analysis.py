# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "numpy",
#     "matplotlib",
# ]
# ///

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import glob
import os

def _setup_plotting():
    """
    Loads all modules for plotting. NOTE: Do not use on the HPC Cluster.
    """
    import seaborn as sns
    import matplotlib as mpl
    mpl.use('TkAgg')
    import scienceplots
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.style.use(['science','nature'])
    return plt

class Lifetimes:
    """
    Defines a Class for analysing a System
    based on the MB Lifetimes.
    """
    def __init__(self,
                 path:str="*.dat"):
        """
        Iniates the Class by reading every ".dat" Input File in the given 
        Path.

        Parameters
        ----------
        path : str
            Path to the Input File.
        """
        self.inputs = [f for f in glob.glob(path) if f.endswith(".dat")]
        self.mbs = [MBAnalysis(path=f) for f in self.inputs]
        self.eis = []
        self.tau_mean = []
        self.tau_std = []
        self.dist = []
        self.m_tau = 0.0
        return
    
    @classmethod
    def read_file(cls,
                  path:str) -> "Lifetimes":
        """
        Initiates the Class by reading a given Result File for Plotting.

        Parameters
        ----------
        path : str
            Path to the Input File.
        """
        data = np.loadtxt(path)
        life = cls.__new__(cls)
        life.inputs = []
        life.mbs = []

        life.eis = data[:,0].tolist()
        life.tau_mean = data[:,1].tolist()
        life.tau_std = data[:,2].tolist()
        life.dist = data[:,3].astype(int).tolist()

        life.m_tau = np.sum(np.array(life.tau_mean)*np.array(life.dist)) / np.sum(life.dist)
        return life
    
    @staticmethod
    def _read_means(path:str) -> list:
        """
        Plots the mean tau for given temperatures from the path. 
        NOTE: Needs files called "Lifetimes.dat".

        Parameters
        ----------
        path : str
            Path to the Input Files.
        """
        import re
        from collections import defaultdict

        files = glob.glob(path)
        temp_tau_dict = defaultdict(list)
        
        def extract_temperature(filepath:str) -> float:
            match = re.search(r'T(\d+\.\d+)',filepath)
            if match:
                return float(match.group(1))
            else:
                raise ValueError
            
        def extract_tau(filepath:str) -> float:
            with open(filepath,"r") as f:
                first_line = f.readline().strip()
            match = re.search(r'<tau>_global\s*=\s*([\d.e+\-]+)',first_line)
            if match:
                return float(match.group(1))
            else:
                raise ValueError
            
        for file in files:
            temp = extract_temperature(file)
            tau = extract_tau(file)
            temp_tau_dict[temp].append(tau)

        temp_list = []
        tau_mean_list = []
        tau_std_list = []

        for temperatures in sorted(temp_tau_dict.keys()):
            tau_values = np.array(temp_tau_dict[temperatures])
            temp_list.append(temperatures)
            tau_mean_list.append(np.mean(tau_values))
            tau_std_list.append(np.std(tau_values))

        temp = np.array(temp_list)
        tau_mean = np.array(tau_mean_list)
        tau_std = np.array(tau_std_list)

        return temp,tau_mean,tau_std

    @staticmethod
    def plot_taus(path:str):
        """
        Plots the mean value of Tau vs the temperature.

        Parameters
        ----------
        path : str
            Path to Input Files.
        """
        plt = _setup_plotting()
        temp, tau_mean, tau_std = Lifetimes._read_means(path)

        fig = plt.figure(figsize=(3.39,2.0))
        ax = fig.add_subplot(121)

        ax.errorbar(temp,tau_mean,yerr=tau_std,
                     marker="o",markersize=4,linestyle="None",
                     capsize=3,capthick=1)
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.set_ylabel(r"$\langle \tau (e_{IS},T)\rangle$")
        ax.set_xlabel(r"$T$")

        ax = fig.add_subplot(122)

        tau_std_log = tau_std / tau_mean

        ax.errorbar(1/temp,np.log(tau_mean),yerr=tau_std_log,
                    marker="o",markersize=4,linestyle="None",
                    capsize=3,capthick=1)
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.set_ylabel(r"$\ln \langle \tau (e_{IS},T)\rangle$")
        ax.set_xlabel(r"$1/T$")

        plt.tight_layout()
        plt.savefig("dist.png")
        plt.savefig("dist.pdf")
        plt.show()
        plt.close()
        return
    
    @staticmethod
    def calc_tau_vs_temp(path:str, energies:list=None,thresh:float=1.0,
                         mode:str="Fit") -> dict:
        """
        Calculates the mean tau from Simulations at different Temperatures 
        for given Energies by fitting the underlying Data.

        Parameters
        ==========
        path : str
            Path to Input Files (expanded via globbing).
        energies : list
            List of Energies to investigate.
        thresh : float
            Threshold Value for creating a Cutoff for Fitting.
        """
        import re
        from collections import defaultdict

        files = glob.glob(path)
        
        def extract_temperature(filepath:str) -> float:
            match = re.search(r'T(\d+\.\d+)',filepath)
            if match:
                return float(match.group(1))
            else:
                raise ValueError
            
        if energies is None:
            for file in files:
                try:
                    life_first = Lifetimes.read_file(files[0])
                    if life_first.filter_tau(thresh=thresh):
                        e_min, e_max = min(life_first.eis),max(life_first.eis)
                        energies = np.linspace(e_min,e_max,10)
                        break
                except Exception as e:
                    print(f"Skipping {file}: {e}")
                    continue
        
        if mode == "Fit":
            tau_data = defaultdict(lambda: {
                "temps": [],
                "tau_ln": [], "tau_ln_err": [],
                "tau": [], "tau_err": []
            })
        elif mode == "Bin":
            tau_data = defaultdict(lambda: {
                "temps": [],
                "tau": [], "tau_std": []
            })

        for file in files:
            temp = extract_temperature(file)
            life = Lifetimes.read_file(file)

            if not life.filter_tau(thresh=thresh):
                print(f"    Skipping T={temp:.4f}: {file}")
                print(f"    Filtering failed.")
                continue

            if mode == "Fit":
                life.fit_ln_tau()
                life.fit_tau()

                e_eval, tau_ln, tau_ln_err, tau, tau_err = life.calc_tau(
                    energies=energies
                )

                for i,e in enumerate(e_eval):
                    tau_data[e]["temps"].append(temp)

                    if tau_ln is not None:
                        tau_data[e]["tau_ln"].append(tau_ln[i])
                        tau_data[e]["tau_ln_err"].append(tau_ln_err[i])
                    else:
                        tau_data[e]["tau_ln"].append(np.nan)
                        tau_data[e]["tau_ln_err"].append(np.nan)

                    if tau is not None:
                        tau_data[e]["tau"].append(tau[i])
                        tau_data[e]["tau_err"].append(tau_err[i])
                    else:
                        tau_data[e]["tau"].append(np.nan)
                        tau_data[e]["tau_err"].append(np.nan)

            elif mode == "Bin":
                e_eval, tau, tau_std, counts = life.binned_tau(
                    energies=energies
                )

                for i,e in enumerate(e_eval):
                    tau_data[e]["temps"].append(temp)

                    if tau is not None:
                        tau_data[e]["tau"].append(tau[i])
                        tau_data[e]["tau_std"].append(tau_std[i])
                    else:
                        tau_data[e]["tau"].append(np.nan)
                        tau_data[e]["tau_std"].append(np.nan)

        return tau_data
        
    @staticmethod
    def plot_tau_vs_temp(tau_data:dict,mode:str="Fit"):
        """
        Plots the mean tau vs the Temperature for a given Data-Frame 
        created in calc_tau_vs_temp().
        
        Parameters
        ==========
        tau_data : dict
            Data from calc_tau_vs_temp()
        """
        plt = _setup_plotting()
        fig = plt.figure(figsize=(3.39,2.0))
        ax = fig.add_subplot(121)

        energies_sorted = sorted(tau_data.keys())
        colors = plt.cm.viridis(np.linspace(0,1,len(energies_sorted)))

        for i,e in enumerate(energies_sorted):
            data = tau_data[e]

            sorted_indices = np.argsort(data["temps"])
            temps = np.array(data["temps"])[sorted_indices]

            if mode == "Fit":
                taus = np.array(data["tau_ln"])[sorted_indices]
                errs = np.array(data["tau_ln_err"])[sorted_indices]

            elif mode == "Bin":
                taus = np.array(data["tau"])[sorted_indices]
                errs = np.array(data["tau_std"])[sorted_indices]

            valid = ~np.isnan(taus)

            ax.errorbar(1/temps[valid],taus[valid],yerr=errs[valid],
                        marker="o",markersize=4,linestyle="--",
                        capsize=3,capthick=1,
                        label=f"$e_{{IS}}={e:.3f}$",
                        color=colors[i])
            
        ax.set_yscale("log")
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$\langle \tau (e_{IS},T) \rangle$")

        ax = fig.add_subplot(122)

        for i,e in enumerate(energies_sorted):
            data = tau_data[e]

            sorted_indices = np.argsort(data["temps"])
            temps = np.array(data["temps"])[sorted_indices]

            if mode == "Fit":
                taus = np.array(data["tau"])[sorted_indices]
                errs = np.array(data["tau_err"])[sorted_indices]

            elif mode == "Bin":
                taus = np.array(data["tau"])[sorted_indices]
                errs = np.array(data["tau_std"])[sorted_indices]

            valid = ~np.isnan(taus)

            ax.errorbar(1/temps[valid],taus[valid],yerr=errs[valid],
                        marker="o",markersize=4,linestyle="--",
                        capsize=3,capthick=1,
                        label=f"$e_{{IS}}={e:.3f}$",
                        color=colors[i])
        
        ax.set_yscale("log")
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$\langle \tau (e_{IS},T) \rangle$")
            
        plt.tight_layout()
        plt.savefig("dist.png")
        plt.savefig("dist.pdf")
        plt.show()
        plt.close()
        return
    
    @staticmethod
    def write_tau_vs_temp(tau_data:dict,path:str,mode:str="Fit"):
        """
        Writes the mean tau vs the Temperature for a given Data-Frame 
        created in calc_tau_vs_temp().
        
        Parameters
        ==========
        tau_data : dict
            Data from calc_tau_vs_temp()
        """
        with open(path,"w") as f:
            f.write(f"# <tau> vs. T\n")
            if mode == "Fit":
                f.write(f"# e_IS T tau_ln tau_ln_err tau tau_err\n")
            elif mode == "Bin":
                f.write(f"# e_IS T tau tau_err\n")
            for e in sorted(tau_data.keys()):
                data = tau_data[e]
                sorted_indices = np.argsort(data["temps"])
                for idx in sorted_indices:
                    temp = data["temps"][idx]
                    if mode == "Fit":
                        tau_ln = data["tau_ln"][idx]
                        tau_ln_err = data["tau_ln_err"][idx]
                        tau = data["tau"][idx]
                        tau_err = data["tau_err"][idx]
                        f.write(f"{e:.10f} {temp:.10f} {tau_ln:.10e} {tau_ln_err:.10e} {tau:.10e} {tau_err:.10e}\n")
                    elif mode == "Bin":
                        tau = data["tau"][idx]
                        tau_err = data["tau_std"][idx]
                        f.write(f"{e:.10f} {temp:.10f} {tau:.10e} {tau_err:.10e}\n")
        return
    
    def calc_lifetimes(self,
                       author:str="Rehwald"):
        """
        Performs the MB Algorithm according to the given author.

        Parameters
        ----------
        author : str
            Which Author's MB Algorithm will be used 
            (Doliwa (1 and 2), Rehwald, Buechner).
        """
        for mb in self.mbs:
            mb.mb_rehwald()
        return
    
    def calc_dist(self):
        """
        Calculates the Distribution of MB energies.
        """
        eis_list = []
        tau_list = []
        dist = []
        for mb in self.mbs:
            for i,eis in enumerate(mb.uniq):
                tau = mb.last[i] - mb.first[i]
                if eis in eis_list:
                    idx = eis_list.index(eis)
                    tau_list[idx].append(tau)
                    dist[idx] += 1
                else:
                    eis_list.append(eis)
                    tau_list.append([tau])
                    dist.append(1)
        self.eis = eis_list
        self.dist = dist
        self.tau_mean = [float(np.mean(tlist)) for tlist in tau_list]
        self.tau_std = [float(np.std(tlist)) for tlist in tau_list]
        return
    
    def filter_tau(self,thresh:float=1.0,sim:float=None) -> bool:
        """
        Filters the underlying tau Data by applying a Cutoff
        calculated from the given Threshold and Simulation-Time.

        Parameters
        ==========
        thresh : float
            Threshold Value for calculating the Cutoff.
        sim : float
            Simulation Time in MD Steps.
        """
        if not self.eis:
            print("ERROR: No data to filter. Run calc_dist() first.")
            return False
        
        if not self.tau_mean:
            print("ERROR: tau_mean is empty.")
            return False
        
        try:
            if sim is None:
                cutoff_time = thresh*max(self.tau_mean)
            else:
                cutoff_time = thresh*sim
        except ValueError:
            print("ERROR: Cannot calculate cutoff - tau_mean contains no valid values.")
            return False
        
        valid_mb_data = []

        for i,tau in enumerate(self.tau_mean):
            if 0 < tau <= cutoff_time:
                valid_mb_data.append((self.eis[i],tau))

        if not valid_mb_data:
            print(f"WARNING: No valid MBs found with 0 < tau <= {cutoff_time:.2e}")
            print(f"    tau_mean range: [{min(self.tau_mean):.2e}, {max(self.tau_mean):.2e}]")
            print(f"    All {len(self.tau_mean)} MBs are outside filter range.")
            print(f"    Skipping filter for this dataset.")
            return False
        
        max_tau = max(tau for eis,tau in valid_mb_data)
        e_max_tau = [eis for eis,tau in valid_mb_data if tau == max_tau]
        e_cutoff = min(e_max_tau)

        filtered_eis = []
        filtered_tau_mean = []
        filtered_tau_std = []
        filtered_dist = []

        for i,eis in enumerate(self.eis):
            if eis >= e_cutoff and self.tau_mean[i] > 0:
                filtered_eis.append(eis)
                filtered_tau_mean.append(self.tau_mean[i])
                filtered_tau_std.append(self.tau_std[i])
                filtered_dist.append(self.dist[i])

        self.eis = filtered_eis
        self.tau_mean = filtered_tau_mean
        self.tau_std = filtered_tau_std
        self.dist = filtered_dist

        self.calc_mean()
        return True
    
    def binned_tau(self,
                   energies:list):
        """
        
        """
        if energies is None:
            e_eval = np.linspace(min(self.eis),max(self.eis),10)
        else:
            e_eval = np.array(energies)
        
        n_bins = len(e_eval)-1

        self.e_eval = e_eval

        if n_bins < 1:
            print("ERROR: Need at least two energy values to create bins.")
            return None, None, None, None
        
        bin_taus = [[] for _ in range(n_bins)]

        for i,eis in enumerate(self.eis):
            tau = self.tau_mean[i]
            dist = self.dist[i]

            for b in range(n_bins):
                if e_eval[b] <= eis < e_eval[b+1]:
                    bin_taus[b].extend([tau] * dist)
                    break
            else:
                if eis == e_eval[-1]:
                    bin_taus[-1].extend([tau] * dist)
        
        bin_centers = []
        tau_mean_bin = []
        tau_std_bin = []
        counts = []

        for b in range(n_bins):
            bin_center = (e_eval[b]+e_eval[b+1])/2
            bin_centers.append(bin_center)

            if len(bin_taus[b]) > 0:
                tau_mean_bin.append(np.mean(bin_taus[b]))
                tau_std_bin.append(np.std(bin_taus[b]))
                counts.append(len(bin_taus[b]))
            else:
                tau_mean_bin.append(np.nan)
                tau_std_bin.append(np.nan)
                counts.append(0)
        
        self.binned_energies = np.array(bin_centers)
        self.binned_tau_mean = np.array(tau_mean_bin)
        self.binned_tau_std = np.array(tau_std_bin)
        self.binned_counts = np.array(counts)

        return (self.binned_energies, self.binned_tau_mean, self.binned_tau_std,
                self.binned_counts) 
    
    def fit_tau(self):
        """
        Fits the tau(e_IS) values as an exponential, quadratic function:
        exp(a+b*e*c*e**2)
        """
        from scipy.optimize import curve_fit
        
        def exp_quadratic(e,a,b,c):
            exponent = a+b*e+c*e**2
            exponent = np.clip(exponent,-100,100)
            if np.any(np.isinf(np.exp(exponent))):
                print("WARNING: Overflow in exponential fit!")
            return np.exp(exponent)

        if not hasattr(self, "a_ln_tau"):
            self.fit_ln_tau()

        p0 = [self.a_ln_tau,self.b_ln_tau,self.c_ln_tau]

        dist_array = np.array(self.dist,dtype=float)
        sigma = np.array(self.tau_mean) / np.sqrt(dist_array)

        try:
            popt, pcov = curve_fit(
                exp_quadratic,
                self.eis,
                self.tau_mean,
                p0=p0,
                sigma=sigma,
                absolute_sigma=False,
                maxfev=10000000
            )
            self.a_tau = popt[0]
            self.b_tau = popt[1]
            self.c_tau = popt[2]
            self.fit_tau_cov = pcov

            return popt, pcov
        except Exception as e:
            print(f"ERROR: Fitting failed with: {type(e).__name__}: {e}")
            return None,None

    def fit_ln_tau(self):
        """
        Fits the ln[tau(e_IS)] values as an quadratic function:
        a+b*e*c*e**2
        """
        from scipy.optimize import curve_fit
        
        def linear_quadratic(e,a,b,c):
            return (a+b*e+c*e**2)
        
        # Since there are MBs occuring for only one timestep (tau=0.0), 
        # we need to filter them.

        tau_array = np.array(self.tau_mean)
        valid_mask = tau_array > 0

        eis_valid = np.array(self.eis)[valid_mask]
        tau_valid = tau_array[valid_mask]
        dist_valid = np.array(self.dist)[valid_mask]

        p0 = [0.0,1.0,0.0]

        ln_tau_mean = np.log(tau_valid)
        sigma = 1.0 / np.sqrt(dist_valid)

        try:
            popt, pcov = curve_fit(
                linear_quadratic,
                eis_valid,
                ln_tau_mean,
                p0=p0,
                sigma=sigma,
                absolute_sigma=False,
                maxfev=50000
            )
            self.a_ln_tau = popt[0]
            self.b_ln_tau = popt[1]
            self.c_ln_tau = popt[2]
            self.fit_ln_tau_cov = pcov
            return popt,pcov
        except Exception as e:
            print(f"ERROR: Fitting failed with: {type(e).__name__}: {e}")
            return None,None  

    def calc_tau(self,energies:list=None):
        """
        Calculates the mean tau from Simulations
        for given Energies by fitting the underlying Data.

        Parameters
        ==========
        energies : list
            List of Energies to investigate.
        """
        if not hasattr(self,"a_ln_tau") and not hasattr(self,"a_tau"):
            self.fit_ln_tau()
            self.fit_tau()

        if energies is None:
            e_eval = np.linspace(min(self.eis),max(self.eis),10)
        else:
            e_eval = np.array(energies)

        cov_ln = self.fit_ln_tau_cov
        tau_ln = np.exp(self.a_ln_tau + self.b_ln_tau*e_eval + self.c_ln_tau*e_eval**2)
        tau_ln_err = np.zeros_like(e_eval)

        for i, e in enumerate(e_eval):
            jacobian = np.array([1,e,e**2])
            var_linear = jacobian @ cov_ln @ jacobian
            tau_ln_err[i] = tau_ln[i] * np.sqrt(var_linear)

        cov = self.fit_tau_cov
        tau = np.exp(self.a_tau + self.b_tau*e_eval + self.c_tau*e_eval**2)
        tau_err = np.zeros_like(e_eval)

        for i, e in enumerate(e_eval):
            jacobian = np.array([1,e,e**2])
            var_linear = jacobian @ cov @ jacobian
            tau_err[i] = tau[i] * np.sqrt(var_linear)

        return e_eval, tau_ln, tau_ln_err, tau, tau_err
    
    def write_fits(self,
                   path:str,
                   energies:list=None):
        """
        Writes the mean tau from Simulations
        for given Energies gained by fitting the underlying Data.

        Parameters
        ==========
        path : str
            Path to Output File.
        energies : list
            List of Energies to investigate.
        """
        e_eval,tau_ln,tau_ln_err,tau,tau_err = self.calc_tau(energies=energies)
        with open(path,"w") as f1:
            f1.write(f"# ln fit: a = {self.a_ln_tau:.10e}, b = {self.b_ln_tau:.10e}, c = {self.c_ln_tau:.10e}\n")
            f1.write(f"# fit: a = {self.a_tau:.10e}, b = {self.b_tau:.10e}, c = {self.c_tau:.10e}\n")
            f1.write(f"# e_IS tau_ln tau_ln_err tau tau_err\n")
            for i,e in enumerate(e_eval):
                f1.write(f"{e:.10f} {tau_ln[i]:.10e} {tau_ln_err[i]:.10e} {tau[i]:.10e} {tau_err[i]:.10e}\n")
        return

    def calc_mean(self):
        """
        Calculates the Mean of Waiting Time from the MB Distribution.
        """
        m_tau = 0.0
        for i,eis in enumerate(self.eis):
            m_tau += self.tau_mean[i]*self.dist[i]
        self.m_tau = m_tau/len(self.eis)
        return
    
    def plot_dist(self):
        """
        Plots the Distribution of MBs and there Waiting Times.
        """
        plt = _setup_plotting()
        fig = plt.figure(figsize=(3.39,2.0))
        ax = fig.add_subplot(121)

        p1 = ax.plot(self.eis,self.tau_mean,marker="o",markersize=4,linestyle="None")
        ax.axhline(self.m_tau,linestyle="--")

        if hasattr(self,"a_tau"):
            e_fit = np.linspace(min(self.eis),max(self.eis),200)
            tau_fit = np.exp(self.a_tau+self.b_tau*e_fit+self.c_tau*e_fit**2)
            ax.plot(e_fit,tau_fit,"--",
                    label=r"$\exp(a+be_{IS}+ce_{IS}^2)$")
            
        if hasattr(self,"a_ln_tau"):
            e_fit = np.linspace(min(self.eis),max(self.eis),200)
            tau_fit = np.exp(self.a_ln_tau+self.b_ln_tau*e_fit+self.c_ln_tau*e_fit**2)
            ax.plot(e_fit,tau_fit,"--",
                    label=r"$a+be_{IS}+ce_{IS}^2$")
            
        if hasattr(self,"binned_energies"):
            ax.plot(self.binned_energies,self.binned_tau_mean,
                    marker="x",markersize=5,linestyle="None")
            for i in self.e_eval:
                ax.vlines(i,min(self.tau_mean),max(self.tau_mean))

        ax.set_yscale("log")
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.set_xlabel(r"$e_{IS}$")
        ax.set_ylabel(r"$\langle \tau (e_{IS})\rangle$")

        ax.set_ylim(min(self.tau_mean),max(self.tau_mean))

        ax.legend()

        ax = fig.add_subplot(122)

        p2 = ax.plot(self.eis,self.dist,marker="o",markersize=4,linestyle="None")

        ax.set_xlabel(r"$e_{IS}$")
        ax.set_ylabel(r"$\varphi (e_{IS})$")
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        plt.tight_layout()
        plt.savefig("dist.png")
        plt.savefig("dist.pdf")
        plt.show()
        plt.close()
        return
    
    def write_results(self,
                      path:str) -> None:
        """
        Writes the MB lifetime statistics to an Output File.

        Parameters
        ----------
        path : str
            Path to Output File.
        """
        with open(path,"w") as f1:
            f1.write(f"# <tau>_global = {self.m_tau:.8e}\n")
            for e,tau_m,tau_s,d in zip(
                self.eis,self.tau_mean,self.tau_std,self.dist
            ):
                f1.write(f"{e:.10f} {tau_m:.8e} {tau_s:.8e} {d:d}\n")
        return
    
    def write_binned(self,
                     path:str,
                     energies:list=None) -> None:
        """
        
        """
        binned_energies,tau_mean,tau_std,counts = self.binned_tau(energies=energies)
        with open(path,"w") as f1:
            f1.write(f"# e_IS tau_mean tau_std counts\n")
            for i,e in enumerate(binned_energies):
                f1.write(f"{e:.10f} {tau_mean[i]:.10e} {tau_std[i]:.10e} {counts[i]:.10e}\n")
        return
        

class IBMAnalysis:
    """
    Defines a Class for using the Interval Bisection Method 
    on a LAMMPS Trajectory.
    """
    def __init__(self,
                 path:str):
        """
        Initiates the Class by reading a given Input File.

        Parameters
        ----------
        path : str
            Path to the Input File.
        """
        self.times,self.e_is = (np.atleast_1d(x) for x in np.loadtxt(path,unpack=True))
        sort = np.argsort(self.times)
        self.times = self.times[sort]
        self.e_is = self.e_is[sort]
        self.e_is = np.round(self.e_is,16)
        self.path = path

        _, unique_idx = np.unique(self.times[::-1],return_index=True)
        unique_idx = len(self.times) - 1 - unique_idx
        unique_idx.sort()
        self.times = self.times[unique_idx]
        self.e_is = self.e_is[unique_idx]
        return
    
    def calc_jumps(self) -> None:
        """
        Prints the Timepoints in between the Energy jumps.
        """
        for i, eis in enumerate(self.e_is):
            if i == 0:
                continue
            if (self.e_is[i-1] != eis) and (abs(self.times[i-1]-self.times[i]) != 1):
                print(f"{int(self.times[i-1])} {int(self.times[i])}")
        return
    
    def clean_file(self) -> None:
        """
        Cleans the given Input File from unnecessary Timepoints.
        """
        keep = set()
        keep.add(0)
        keep.add(len(self.e_is)-1)
        for i,eis in enumerate(self.e_is):
            if i == 0 or i == (len(self.e_is)-1):
                continue
            if (self.e_is[i-1] != eis):
                keep.add(i)
                keep.add(i-1)
        
        keep = sorted(keep)

        times_clean = self.times[keep]
        eis_clean = self.e_is[keep]
        with open(self.path, "w") as f:
            for t,e in zip(times_clean, eis_clean):
                f.write(f"{int(t)} {e:.15f}\n")
        return

class MBAnalysis:
    """
    Defines a Class for using the Interval Bisection Method 
    on a LAMMPS Trajectory.
    """
    def __init__(self,
                 path:str):
        """
        Initiates the Class by reading a given Input File.

        Parameters
        ----------
        path : str
            Path to the Input File.
        """
        self.times,self.e_is = (np.atleast_1d(x) for x in np.loadtxt(path,unpack=True))
        sort = np.argsort(self.times)
        self.times = self.times[sort]
        self.e_is = self.e_is[sort]
        self.e_is = np.round(self.e_is,16)
        self.first = None
        self.last = None
        self.uniq = None
        self.doliwa = None
        self.nr_int = 0
        return
    
    @classmethod
    def read_file(cls,
                  path:str):
        """
        Initiates the Class by reading an already written MB File.

        Parameters
        ----------
        path : str
            Path to Input File.
        """
        data = np.loadtxt(path)
        
        if data.ndim == 1:
            data = data[None, :]
        
        mb = cls.__new__(cls)
        mb.times = None
        mb.e_is = None

        mb.first = data[:,0]
        mb.last = data[:,1]
        mb.uniq = data[:,2]
        mb.nr_int = len(mb.uniq)
        return mb
    
    def mb_doliwa(self,
                  version:int=0) -> None:
        """
        Gathers the Functions necessary for performing the MB Algorithm
        according to B. Doliwa (2002).

        Parameters
        ----------
        version : int
            Defines the choosen version of Doliwa's algorithm (0 or 1).
        """
        self.calc_intervals()
        self.calc_cut_doliwa(version=version)
        self.calc_combine()
        self.calc_remove()
        self._check_duplicates()
        self._check_overlaps()
        return

    def mb_buechner(self,
                    alpha:float=0.1) -> None:
        """
        Gathers the Functions necessary for performing the MB Algorithm
        according to S. Büchner (2002).

        Parameters
        ----------
        alpha : float
            The maximal Time a Interval needs to be,
            to be classified as a single Basin. (Usually tau/10.)
        """
        self.calc_intervals()
        self.calc_combine_buechner(alpha=alpha)
        self.calc_remove()
        self._check_duplicates()
        self._check_overlaps()
        return

    def mb_rehwald(self) -> None:
        """
        Gathers the Functions necessary for performing the MB Algorithm
        according to C. Rehwald (2012).
        """
        self.calc_intervals()
        self.calc_cut_rehwald()
        self.calc_cut_rehwald()
        self.calc_combine()
        self.calc_remove()
        self._check_duplicates()
        self._check_overlaps()
        return

    def calc_intervals(self) -> None:
        """
        Calculates the Intervals given by the first and last Occurance
        of a Minimum of the Trajectory and stores them in uniq, first and last. 
        """
        uniq = np.unique(self.e_is)
        first_list = []
        last_list = []
        uniq_list = []
        for eis in uniq:
            idx = np.where(self.e_is == eis)[0]
            first = self.times[idx[0]]
            last = self.times[idx[-1]]
            first_list.append(first)
            last_list.append(last)
            uniq_list.append(eis)
        self.first = np.array(first_list)
        self.last = np.array(last_list)
        self.uniq = np.array(uniq_list)
        self.nr_int = len(uniq_list)
        self._sort_intervals()
        return
    
    def calc_cut_rehwald(self) -> None:
        """
        Calculates the Cut of Intervals overlapping less than 50%
        according to C. Rehwald (2012).
        """
        updates = {}
        for i in range(self.nr_int):
            for j in range(i+1,self.nr_int):
                if self.first[i] < self.first[j] < self.last[i] < self.last[j]:
                    overlap = self._calc_overlap(i,j)
                    if overlap <= 0.5:
                        updates[j] = max(updates.get(j,self.first[j]),self.last[i])
        if len(updates) > 0:
            for j, new_first in updates.items():
                self.first[j] = new_first
        self._sort_intervals()
        return

    def calc_cut_doliwa(self,
                        version:int=0) -> None:
        """
        Calculates the Cut of Intervals overlapping less than 50%
        according to B. Doliwa (2002) in two different interpretations.

        Parameters
        ----------
        version : int
            Version of the cutting Algorithm.
        """
        """
        updates_first = {}
        updates_last = {}
        for i in range(self.nr_int):
            for j in range(i+1,self.nr_int):
                if self.first[i] < self.first[j] < self.last[i] < self.last[j]:
                    overlap = self.calc_overlap(i,j)
                    if overlap <= 0.5:
                        roll = np.random.randint(0,2)
                        if roll == 0:
                            if version == 0:
                                updates_last[i] = ("traj",j)
                            else:
                                updates_last[i] = min(updates_last.get(i,self.last[i]),self.first[j])
                        else:
                            if version == 0:
                                updates_first[j] = ("traj",i)
                            else:
                                updates_first[j] = max(updates_first.get(j,self.first[j]),self.last[i])
        for i,val in updates_last.items():
            if isinstance(val,tuple):
                _,j = val
                self.cut_doliwa_traj(i=i,j=j)
            else:
                self.last[i] = val
        for j,val in updates_first.items():
            if isinstance(val, tuple):
                _,i = val
                self.cut_doliwa_traj(i=j,j=i)
            else:
                self.first[j] = val
        self.sort_intervals()
        return
        """
        pass
    
    def cut_doliwa_traj(self,
                        i:int,
                        j:int) -> None:
        """
        Performs the first Version of Cutting by B. Doliwa:
        Cuts Interval (i) to the first or last Occurance in the 
        original Trajectory outside of the second Interval (j).

        Parameters
        ----------
        i : int
            Index of the Interval to cut.
        j : int
            Second Interval.
        """
        """
        eps = self.doliwa[i]
        if self.first[i] < self.first[j]:
            mask = (self.times <= self.first[j]) & (self.e_is == eps)
            idx = np.where(mask)[0]
            if idx.size > 0:
                self.last[i] = self.times[idx.min()]
        else:
            mask = (self.times >= self.last[j]) & (self.e_is == eps)
            idx = np.where(mask)[0]
            if idx.size > 0:
                self.first[i] = self.times[idx.max()]
        return
        """
        pass
    
    def calc_combine(self) -> None:
        """
        Performs the Combination of Intervals with an Overlap of more than 50%.
        """
        parent = list(range(self.nr_int))
        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]
        def union(x,y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py
                return True
            return False
        combine_count = 0
        for i in range(self.nr_int):
            for j in range(i+1,self.nr_int):
                if self.first[i] < self.first[j] < self.last[i] < self.last[j]:
                    overlap = self._calc_overlap(i,j)
                    if overlap > 0.5:
                        if union(i,j):
                            combine_count += 1
        clusters = {}
        for i in range(self.nr_int):
            clusters.setdefault(find(i),[]).append(i)

        new_first = []
        new_last = []
        new_uniq = []
        for cluster_id, members in clusters.items():
            if len(members) == 1:
                idx = members[0]
                new_first.append(self.first[idx])
                new_last.append(self.last[idx])
                new_uniq.append(self.uniq[idx])
            else:
                f = np.min(self.first[members])
                l = np.max(self.last[members])
                min_energy = np.min(self.uniq[members])
                new_first.append(f)
                new_last.append(l)
                new_uniq.append(min_energy)
        self.first = np.array(new_first)
        self.last = np.array(new_last)
        self.uniq = np.array(new_uniq)
        self.nr_int = len(new_first)
        self._sort_intervals()
        return
    
    def calc_combine_buechner(self,
                              alpha:float) -> None:
        """
        Combines Intervals 
        according to S. Büchner's Algorithm (2002).

        Parameters
        ----------
        alpha : float
            The maximal Time a Interval needs to be,
            to be classified as a single Basin. (Usually tau/10.)
        """
        """
        i = 0
        span = []
        while i < len(self.uniq):
            span_first = self.first[i]
            span_last = self.last[i]
            min_energy = self.uniq[i]
            max_energy = self.uniq[i]
            merged = False
            j = i+1
            while j < len(self.uniq):
                if self.first[j] < span_last:
                    span_first = min(span_first,self.first[j])
                    span_last = max(span_last,self.last[j])
                    min_energy = min(min_energy,self.uniq[j])
                    max_energy = max(max_energy,self.uniq[j])
                    self.first = np.delete(self.first,j)
                    self.last = np.delete(self.last,j)
                    self.uniq = np.delete(self.uniq,j)
                    continue
                else:
                    break
                j += 1
            if merged:
                self.first[i] = span_first
                self.last[i] = span_last
                if (span_last - span_first) >= alpha:
                    self.uniq[i] = min_energy
                else:
                    self.uniq[j] = max_energy
            i += 1
        self.sort_intervals()
        return
        """
        pass
    
    def calc_remove(self) -> None:
        """
        Performs the Deletion of Intervals that lie within other Intervals.
        """
        remove = set()
        for i in range(self.nr_int):
            if i in remove:
                continue
            for j in range(self.nr_int):
                if i == j or j in remove:
                    continue
                if self.first[i] <= self.first[j] and self.last[i] >= self.last[j]:
                    self.uniq[i] = min(self.uniq[i],self.uniq[j])
                    remove.add(j)
        keep_mask = np.ones(self.nr_int,dtype=bool)
        keep_mask[list(remove)] = False
        self.first = self.first[keep_mask]
        self.last = self.last[keep_mask]
        self.uniq = self.uniq[keep_mask]
        self.nr_int = np.sum(keep_mask)
        self._sort_intervals()
        return

    def _calc_overlap(self,
                      i:int,
                      j:int) -> float:
        """
        Calculates the Overlap of two Intervals i and j.

        Parameters
        ----------
        i : int
            Index of the first Interval.
        j : int
            Index of the second Interval.
        """
        max_len = max((self.last[i]-self.first[i]),self.last[j]-self.first[j])
        overlap = max(self.first[i],self.first[j])-min(self.last[i],self.last[j])
        return (overlap/max_len)

    def _sort_intervals(self) -> None:
        """
        Sorts the Arrays uniq, first and last in increasing Start-Time.
        """
        sort = np.argsort(self.first)
        self.first = self.first[sort]
        self.last = self.last[sort]
        self.uniq = self.uniq[sort]
        return
    
    def _check_duplicates(self) -> None:
        """
        Checks the intervals for duplicated energies.
        """
        unique_energies = set(self.uniq)
        if len(unique_energies) < self.nr_int:
            counts = {}
            for e in self.uniq:
                counts[e] = counts.get(e,0)+1
            duplicates = {e: count for e, count in counts.items() if count > 1}
            print(f"WARNING: {len(duplicates)} energies appear multiple times.")
        return
    
    def _check_overlaps(self) -> None:
        """
        Checks for overlapping intervals.
        """
        overlaps = []
        for i in range(self.nr_int):
            for j in range(i+1,self.nr_int):
                left = max(self.first[i],self.first[j])
                right = min (self.last[i],self.last[j])
                if right > left:
                    overlaps.append((i,j,right - left))
        if overlaps:
            print(f"WARNING: {len(overlaps)} intervals overlap.")
        return
    
    def plot_eis(self) -> None:
        """
        Plots the initial Trajectory of Energetic Minimas.
        """
        plt = _setup_plotting()
        fig = plt.figure(figsize=(3.39,2.0))
        ax = fig.add_subplot()

        p1 = ax.plot(self.times,self.e_is,marker="o",markersize=4,linestyle="None")

        ax.set_xlabel("Time")
        ax.set_ylabel(r"$e_{IS}$")
        
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        plt.tight_layout()
        plt.savefig("eis.png")
        plt.savefig("eis.pdf")
        plt.show()
        plt.close()
        return

    def plot_intervals(self) -> None:
        """
        Plots the Intervals of Energetic Minimas.
        """
        plt = _setup_plotting()
        fig = plt.figure(figsize=(3.39,2.0))
        ax = fig.add_subplot()

        for eis,start,end in zip(self.uniq,self.first,self.last):
            ax.hlines(y=eis,xmin=start,xmax=end,linewidth=0.8,color="black")

        ax.set_xlabel("Time")
        ax.set_ylabel(r"$e_{IS}$")
        
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        plt.tight_layout()
        plt.savefig("eis.png")
        plt.savefig("eis.pdf")
        plt.show()
        plt.close()
        return
    
    def plot_all(self) -> None:
        """
        Plots the initial IS Trajectory and the resulting MBs.
        """
        plt = _setup_plotting()
        fig = plt.figure(figsize=(3.39,2.0))
        ax = fig.add_subplot()

        p1 = ax.plot(self.times,self.e_is,marker="x",markersize=4,linestyle="None")
        for eis,start,end in zip(self.uniq,self.first,self.last):
            ax.hlines(y=eis,xmin=start,xmax=end,linewidth=0.8,color="black")

        ax.set_xlabel("Time")
        ax.set_ylabel(r"$e_{IS}$")
        
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        plt.tight_layout()
        plt.savefig("eis.png")
        plt.savefig("eis.pdf")
        plt.show()
        plt.close()
        return
    
    def write_results(self,
                      path:str) -> None:
        """
        Writes the resulting Intervals to an Output File.

        Parameters
        ----------
        path : str
            Path to Output File.
        """
        with open(path,"w") as f1:
            for i,eis in enumerate(self.uniq):
                f1.write(f"{self.first[i]} {self.last[i]} {eis}\n")
        return
    
def main():
    """
    :Author: Simon Georg Kellers
    :Email: s_kell14@uni-muenster.de
    :Date: 2026
    :Version: 5.1.0
    
    DESCRIPTION
        Performs different Analyses on minimized IS-Energies:
            IBM : Interval Bisection Method
                Detects Timepoints inbetween the IS-Energy changes.
            MB : Meta-Basin Algorithm
                Identifes meta-Basins based on IS-Energy Intervals,
                Overlap Handling and different Combination Rules.
            Clean
                Cleans the given IS Trajectory File by removing unnecessary
                timepoints.
            Lifetimes
                Calculates the Lifetimes based on the IS-Trajectory.
            TauT
                Calculates the mean Lifetime of given Energies at different
                Temperatures.

    USAGE
        python3 MB_Analysis.py [-h] [-m MODE]

        -h, --help
                prints this Message

        -i, --input INPUT
                Path to Input File.

        -m, --mode MODE
                IBM, MB, Clean, Lifetimes

        -o, --output OUTPUT
                Output for the MB Analysis' Results.

        -a, --author AUTHOR
                which Author's MB Algorithm will be used 
                (Doliwa (1 and 2), Rehwald, Buechner).

        -t, --timepoint TIMEPOINT
                the Timepoint to compare for a Restart File.
    """
    parser = argparse.ArgumentParser(description=(
        "Performs different Analyses on minimized IS-Energies:\n" \
        "   IBM : Interval Bisection Method\n" \
        "        Detects Timepoints inbetween the IS-Energy changes.\n" \
        "    MB : Meta-Basin Algorithm\n" \
        "        Identifes meta-Basins based on IS-Energy Intervals,\n" \
        "        Overlap Handling and different Combination Rules.\n" \
        "    Lifetimes\n" \
        "        Calculates the Lifetimes based on the IS Trajectory.\n" \
        "    Clean\n" \
        "        Shortens the IS Trajectory File and saves only\n" \
        "        necessary Transitions.\n" \
        "    TauT\n" \
        "        Calculates the mean Lifetime of given Energies."
    ),
    epilog="For more Details and Plotting Functions, see the Documentation.",
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-i","--input",type=str,required=True,help="Path to Input File.")
    parser.add_argument("-m","--mode",type=str,
                        choices=["IBM","MB","Clean","Lifetimes","Fitting","TauT",
                                 "MB_Plot","Life_Plot","Tau_Plot","TauT_Plot"],required=True)
    # Parameters for the MB Algorithm
    parser.add_argument("-o","--output",type=str,help="Output for the MB Analysis' Results.")
    parser.add_argument("-a","--author",
                        const="Rehwald",
                        nargs="?",
                        type=str,
                        choices=["Doliwa1","Doliwa2","Rehwald","Buechner"],
                        help="which Author's MB Algorithm will be used.")
    # Parameters for the Analysis
    parser.add_argument("-t","--thresh",type=float,default=1.0,
                        help="Threshold Value for the Fitting of tau.")
    parser.add_argument("-e","--energies",nargs="+",type=float,
                        help="Energies to determine the tau from fits.")
    parser.add_argument("-s","--submode",type=str,default="Fit",
                        help="Mode for the Fitting or Binning of e_IS.")
    
    args = parser.parse_args()

    if args.mode == "IBM":
        mb = IBMAnalysis(path=args.input)
        mb.calc_jumps()

    elif args.mode == "Clean":
        mb = IBMAnalysis(path=args.input)
        mb.clean_file()
    
    elif args.mode == "MB":
        mb = MBAnalysis(path=args.input)
        mb.mb_rehwald()
        mb.write_results(path=args.output)

    elif args.mode == "Lifetimes":
        life = Lifetimes(path=args.input)
        life.calc_lifetimes(author=args.author)
        life.calc_dist()
        life.calc_mean()
        life.write_results(path=args.output)

    elif args.mode == "Fitting":
        life = Lifetimes.read_file(path=args.input)
        life.filter_tau(thresh=args.thresh)
        life.write_fits(path=args.output,energies=args.energies)
    
    elif args.mode == "TauT":
        data = Lifetimes.calc_tau_vs_temp(path=args.input,energies=args.energies,thresh=args.thresh)
        Lifetimes.write_tau_vs_temp(tau_data=data)

    elif args.mode == "MB_Plot":
        mb = MBAnalysis.read_file(path=args.input)
        mb.plot_intervals()

    elif args.mode == "Life_Plot":
        life = Lifetimes.read_file(path=args.input)
        life.filter_tau(thresh=args.thresh)
        life.fit_tau()
        life.fit_ln_tau()
        life.binned_tau(energies=args.energies)
        life.plot_dist()

    elif args.mode == "Tau_Plot":
        Lifetimes.plot_taus(path=args.input)

    elif args.mode == "TauT_Plot":
        data = Lifetimes.calc_tau_vs_temp(path=args.input,
                                          energies=args.energies,
                                          thresh=args.thresh,
                                          mode=args.submode)
        Lifetimes.plot_tau_vs_temp(tau_data=data,
                                   mode=args.submode)

    return

if __name__ == "__main__":
    main()