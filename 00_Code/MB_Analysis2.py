# /// script
# requires-python = ">=3.7"
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
import time
from tqdm import tqdm

# ── FitResult ─────────────────────────────────────────────────────────────────

class FitResult:
    """
    Holds the result of a quadratic fit to ln(tau) vs e_IS:
        ln(tau) = a + b*e + c*e^2
    """
    def __init__(self, a:float, b:float, c:float, cov:float):
        self.a   = a
        self.b   = b
        self.c   = c
        self.cov = cov
        return

    def eval(self, energies:list):
        """
        Evaluates tau and its propagated error at the given energies.

        Parameters
        ----------
        energies : array-like
            Energies at which to evaluate the fit.

        Returns
        -------
        tau : np.ndarray
        tau_err : np.ndarray
        """
        e      = np.asarray(energies)
        ln_tau = self.a + self.b*e + self.c*e**2
        tau    = np.exp(ln_tau)
        errors = np.array([
            tau[i] * np.sqrt(
                np.array([e[i]**2, e[i], 1]) @ self.cov @ np.array([e[i]**2, e[i], 1])
            )
            for i in range(len(e))
        ])
        return tau, errors


# ── Lifetimes ─────────────────────────────────────────────────────────────────

class Lifetimes:
    """
    Holds the MB lifetime statistics for a single simulation.
    """
    def __init__(self, eis:list, tau_mean:list, tau_std:list, dist:list):
        self.eis      = list(eis)
        self.tau_mean = list(tau_mean)
        self.tau_std  = list(tau_std)
        self.dist     = list(dist)
        self.m_tau    = 0.0
        return

    def calc_mean(self) -> None:
        """
        Calculates the global mean waiting time from the MB distribution.
        """
        if not self.eis:
            return
        total_weight = sum(self.dist)
        self.m_tau   = sum(t*d for t,d in zip(self.tau_mean, self.dist)) / total_weight
        return

    def filter_tau(self, thresh:float=1.0, sim:int=None) -> bool:
        """
        Filters MB data by applying a cutoff on tau_mean.

        Parameters
        ----------
        thresh : float
            Multiplier for the cutoff (relative to max tau or sim time).
        sim : float, optional
            Simulation time in MD steps. If None, uses max(tau_mean).

        Returns
        -------
        bool
            False if filtering removes all data or input is invalid.
        """
        if not self.eis or not self.tau_mean:
            print("ERROR: No data to filter.")
            return False

        try:
            cutoff_time = thresh * (sim if sim is not None else max(self.tau_mean))
        except ValueError:
            print("ERROR: Cannot calculate cutoff.")
            return False

        valid = [(self.eis[i], self.tau_mean[i])
                 for i in range(len(self.eis))
                 if 0 < self.tau_mean[i] <= cutoff_time]

        if not valid:
            print(f"WARNING: No valid MBs found with 0 < tau <= {cutoff_time:.2e}")
            return False

        max_tau  = max(tau for _, tau in valid)
        e_cutoff = min(e for e, tau in valid if tau == max_tau)

        idx = [i for i, e in enumerate(self.eis)
               if e >= e_cutoff and self.tau_mean[i] > 0]

        self.eis      = [self.eis[i]      for i in idx]
        self.tau_mean = [self.tau_mean[i] for i in idx]
        self.tau_std  = [self.tau_std[i]  for i in idx]
        self.dist     = [self.dist[i]     for i in idx]

        self.calc_mean()
        return True


# ── LifetimesFitter ───────────────────────────────────────────────────────────

class LifetimesFitter:
    """
    Performs fitting and binning on a Lifetimes object.
    """
    def __init__(self, life:Lifetimes):
        self.life    = life
        self.fit_ln  = None
        self.fit_dir = None

        self.binned_energies = None
        self.binned_tau_mean = None
        self.binned_tau_std  = None
        self.binned_dist     = None

    def fit_ln_tau(self) -> FitResult:
        """
        Fits ln(tau) vs e_IS directly on the raw (unbinned) data,
        weighted by sqrt(dist).

        Returns
        -------
        FitResult or None if the fit fails.
        """
        tau_array  = np.array(self.life.tau_mean)
        valid      = tau_array > 0

        eis_valid  = np.array(self.life.eis)[valid]
        tau_valid  = tau_array[valid]
        dist_valid = np.array(self.life.dist)[valid]

        ln_tau  = np.log(tau_valid)
        weights = np.sqrt(dist_valid)

        try:
            coeffs, cov = np.polyfit(eis_valid, ln_tau, deg=2, w=weights, cov=True)
            self.fit_ln = FitResult(a=coeffs[2], b=coeffs[1], c=coeffs[0], cov=cov)
            return self.fit_ln
        except Exception as e:
            print(f"ERROR: fit_ln_tau failed: {type(e).__name__}: {e}")
            return None

    def fit_dir_tau(self, energies:list=None, nr_bins:int=10,
                    min_count:int=1) -> FitResult:
        """
        Fits ln(tau) vs e_IS on binned data, weighted by sqrt(counts).

        Parameters
        ----------
        energies : list, optional
            Bin edges. If None, nr_bins equally spaced bins are used.
        nr_bins : int
            Number of bins (only used when energies is None).
        min_count : int
            Minimum occupation of a bin to be considered.

        Returns
        -------
        FitResult or None if the fit fails.
        """
        if self.binned_energies is None:
            self.binned_tau(energies=energies, nr_bins=nr_bins,
                            min_count=min_count)

        valid_mask = (
            np.isfinite(self.binned_tau_mean) &
            (self.binned_tau_mean > 0) &
            (self.binned_dist     > 0)
        )
        eis_valid  = self.binned_energies[valid_mask]
        tau_valid  = self.binned_tau_mean[valid_mask]
        dist_valid = self.binned_dist[valid_mask]

        if len(eis_valid) == 0:
            print("ERROR: No valid bins for fit_dir_tau.")
            return None

        ln_tau  = np.log(tau_valid)
        weights = np.sqrt(dist_valid)

        try:
            coeffs, cov = np.polyfit(eis_valid, ln_tau, deg=2, w=weights, cov=True)
            self.fit_dir = FitResult(a=coeffs[2], b=coeffs[1], c=coeffs[0], cov=cov)
            return self.fit_dir
        except Exception as e:
            print(f"ERROR: fit_dir_tau failed: {type(e).__name__}: {e}")
            return None

    def calc_tau(self, energies:list=None, min_count:int=1):
        """
        Evaluates both fits at the given energies.
        Runs the fits first if they have not been run yet.

        Parameters
        ----------
        energies : list, optional
            Energies at which to evaluate. Defaults to 10 points over the data range.
        min_count : int
            Minimum occupation for abin to be considered.

        Returns
        -------
        e_eval, tau_ln, tau_ln_err, tau_dir, tau_dir_err
        """
        if self.binned_energies is None:
            self.binned_tau(energies=energies,min_count=min_count)

        if self.fit_ln  is None:
            self.fit_ln_tau()
        if self.fit_dir is None:
            self.fit_dir_tau(energies=energies,min_count=min_count)

        e_eval = self.binned_energies

        tau_ln,  tau_ln_err  = self.fit_ln.eval(e_eval)  if self.fit_ln  is not None else (None, None)
        tau_dir, tau_dir_err = self.fit_dir.eval(e_eval) if self.fit_dir is not None else (None, None)

        return e_eval, tau_ln, tau_ln_err, tau_dir, tau_dir_err

    def binned_tau(self, energies:list=None, nr_bins:int=100,
                   min_count:int=1):
        """
        Calculates the mean tau within bins defined by the given energies.

        Parameters
        ----------
        energies : list, optional
            Bin edges. If None, nr_bins equally spaced bins are used.
        nr_bins : int
            Number of bins (only used when energies is None).
        min_count : int
            Minimum occupation for a bin to be considered.

        Returns
        -------
        bin_centers, tau_mean_bin, tau_std_bin, counts
        """
        e_eval = (np.linspace(min(self.life.eis), max(self.life.eis), nr_bins+1)
                  if energies is None else np.asarray(energies))

        n_bins = len(e_eval) - 1
        if n_bins < 1:
            print("ERROR: Need at least two energy values to define bins.")
            return None, None, None, None

        bin_taus = [[] for _ in range(n_bins)]

        for i, eis in enumerate(self.life.eis):
            tau  = self.life.tau_mean[i]
            dist = self.life.dist[i]
            for b in range(n_bins):
                if e_eval[b] <= eis < e_eval[b+1]:
                    bin_taus[b].extend([tau] * dist)
                    break
            else:
                if eis == e_eval[-1]:
                    bin_taus[-1].extend([tau] * dist)

        bin_centers  = []
        bin_tau_mean = []
        bin_tau_std  = []
        bin_dist     = []

        for b in range(n_bins):
            bin_centers.append((e_eval[b] + e_eval[b+1]) / 2)
            if bin_taus[b]:
                count = len(bin_taus[b])
                bin_dist.append(len(bin_taus[b]))
                if count >= min_count:
                    bin_tau_mean.append(np.mean(bin_taus[b]))
                    bin_tau_std.append(np.std(bin_taus[b]))
                else:
                    bin_tau_mean.append(np.nan)
                    bin_tau_std.append(np.nan)
            else:
                bin_tau_mean.append(np.nan)
                bin_tau_std.append(np.nan)
                bin_dist.append(0)

        self.binned_energies = np.array(bin_centers)
        self.binned_tau_mean = np.array(bin_tau_mean)
        self.binned_tau_std  = np.array(bin_tau_std)
        self.binned_dist     = np.array(bin_dist)

        return (self.binned_energies, self.binned_tau_mean,
                self.binned_tau_std,  self.binned_dist)
    
    def check_binned_tau(self,nr_bins:int=100,min_count:int=1,
                         energies:list=None):
        """
        
        """
        if energies is not None:
            edges = np.asarray(energies)
            n_bins = len(edges)-1
        else:
            e_min = min(self.life.eis)
            e_max = max(self.life.eis)
            edges = np.linspace(e_min,e_max,nr_bins+1)
            n_bins = len(edges)-1

        bin_taus = [[] for _ in range(n_bins)]

        for i,eis in enumerate(self.life.eis):
            tau = self.life.tau_mean[i]
            dist = int(self.life.dist[i])
            for b in range(n_bins):
                if edges[b] <= eis < edges[b+1]:
                    bin_taus[b].extend([tau]*dist)
                    break
            else:
                if eis == edges[-1]:
                    bin_taus[-1].extend([tau]*dist)
        
        centers = (edges[:-1] + edges[1:]) / 2
        counts = np.array([len(bt) for bt in bin_taus])
        tau_mean = np.full(n_bins,np.nan)
        tau_std = np.full(n_bins,np.nan)
        tau_sem = np.full(n_bins,np.nan)
        tau_med = np.full(n_bins,np.nan)
        rel_sem = np.full(n_bins,np.nan)
        skweness = np.full(n_bins,np.nan)

        for b in range(n_bins):
            if counts[b] >= min_count:
                arr = np.array(bin_taus[b])
                m = np.mean(arr)
                s = np.std(arr)
                n = counts[b]

                tau_mean[b] = m
                tau_std[b] = s
                tau_sem[b] = s / np.sqrt(n)
                tau_med[b] = np.median(arr)

                if m > 0:
                    rel_sem[b] = tau_sem[b] / m
                if s > 0 and n >= 3:
                    skweness[b] = np.mean(((arr-m) / s) ** 3)

        return {
            "edges": edges,
            "centers": centers,
            "counts": counts,
            "tau_mean": tau_mean,
            "tau_std": tau_std,
            "tau_sem": tau_sem,
            "tau_med": tau_med,
            "rel_sem": rel_sem,
            "skweness": skweness,
            "bin_taus": bin_taus,
        }
    
    def check_fitted_tau(self,nr_bins:int=50,min_count:int=1,
                         energies:list=None):
        """
        
        """
        stats = self.check_binned_tau(nr_bins=nr_bins,min_count=min_count,
                                      energies=energies)

        centers = stats["centers"]
        tau_m = stats["tau_mean"]
        tau_sem = stats["tau_sem"]
        counts = stats["counts"]

        valid = (np.isfinite(tau_m) & (tau_m > 0)) & (counts >= min_count)

        e_valid = centers[valid]
        tau_valid = tau_m[valid]
        sem_valid = tau_sem[valid]
        cnt_valid = counts[valid]
        n_pts = int(np.sum(valid))

        ln_tau = np.log(tau_valid)

        w_sqrt_n = np.sqrt(cnt_valid)

        sigma_ln = sem_valid / tau_valid
        sigma_ln = np.where(sigma_ln > 0, sigma_ln, np.nanmedian(sigma_ln))

        w_inv_sig = 1.0 / sigma_ln

        w_uniform = np.ones_like(ln_tau)

        results = {}

        for name, weights in [("sqrt_n",w_sqrt_n),
                              ("inv_sig",w_inv_sig),
                              ("uniform",w_uniform)]:
            try:
                coeffs,cov = np.polyfit(e_valid,ln_tau,deg=2,w=weights,cov=True)
            except:
                results[name] = None
                continue

            ln_tau_pred = np.polyval(coeffs,e_valid)
            residuals = ln_tau - ln_tau_pred

            wmean = np.average(ln_tau,weights=weights**2)
            ss_res = np.sum(weights**2 * residuals**2)
            ss_tot = np.sum(weights**2 * (ln_tau-wmean)**2)
            R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

            dof = n_pts - 3
            chi2_red = ss_res / dof if dof > 0 else np.nan

            e_fine = np.linspace(e_valid.min(),e_valid.max(),200)
            ln_tau_fine = np.polyval(coeffs,e_fine)
            tau_fine = np.exp(ln_tau_fine)

            tau_eval = np.exp(ln_tau_pred)
            sigma_fit = np.array([
                tau_eval[i] * np.sqrt(
                    np.array([e_valid[i]**2,e_valid[i],1])
                    @ cov
                    @ np.array([e_valid[i]**2,e_valid[i],1])
                )
                for i in range(n_pts)
            ])

            results[name] = {
                "coeffs":      coeffs,
                "cov":         cov,
                "e_valid":     e_valid,
                "ln_tau":      ln_tau,
                "ln_tau_pred": ln_tau_pred,
                "residuals":   residuals,
                "weights":     weights,
                "sigma_ln":    sigma_ln,
                "counts":      cnt_valid,
                "R2":          R2,
                "chi2_red":    chi2_red,
                "n_pts":       n_pts,
                "e_fine":      e_fine,
                "tau_fine":    tau_fine,
                "tau_eval":    tau_eval,
                "sigma_fit":   sigma_fit,
            }

        return {"bin_stats": stats, "fits": results}
    
    @staticmethod
    def check_activ(tau_data:dict,min_total:int=0) -> dict:
        """
        
        """
        results = {}

        for e, data in tau_data.items():
            sorted_idx = np.argsort(data["temps"])
            temps = np.array(data["temps"])[sorted_idx]
            taus = np.array(data["tau_ln"])[sorted_idx]
            errs = np.array(data["tau_ln_err"])[sorted_idx]
            counts_raw = np.array(
                data.get("counts", [1]*len(data["temps"]))
            )[sorted_idx].astype(float)

            total_counts = int(np.sum(counts_raw))
            if min_total > 0 and total_counts < min_total:
                results[e] = None
                continue

            valid = np.isfinite(taus) & (taus > 0) & (temps > 0)
            if np.sum(valid) < 2:
                results[e] = None
                continue

            inv_T = 1.0 / temps[valid]
            ln_tau = np.log(taus[valid])
            errs_v = errs[valid]
            counts_v = counts_raw[valid]
            n_pts = int(np.sum(valid))

            w_sqrt_n = np.sqrt(counts_v)
            w_uniform = np.ones(n_pts)

            sigma_ln = errs_v / taus[valid]
            median_sig = np.nanmedian(sigma_ln)
            if np.isfinite(median_sig) and median_sig > 0:
                sigma_ln_filled = np.where(sigma_ln > 0, sigma_ln, median_sig)
                w_inv_sig = 1.0 / sigma_ln_filled
            else:
                w_inv_sig = None

            inv_T_fine = np.linspace(inv_T.min(),inv_T.max(),100)

            fit_results = {}

            for name, weights in [("sqrt_n",w_sqrt_n),
                                  ("inv_sig",w_inv_sig),
                                  ("uniform",w_uniform)]:
                try: 
                    coeffs,cov = np.polyfit(inv_T,ln_tau,deg=1,w=weights,cov="unscaled")
                except Exception:
                    fit_results[name] = None
                    continue

                ln_tau_pred = np.polyval(coeffs,inv_T)
                residuals = ln_tau - ln_tau_pred

                wmean = np.average(ln_tau,weights=weights**2)
                ss_res = np.sum(weights**2 * residuals**2)
                ss_tot = np.sum(weights**2 * (ln_tau - wmean)**2)
                R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

                dof = n_pts - 2
                chi2_red = ss_res / dof if dof > 0 else np.nan

                fit_results[name] = {
                    "coeffs": coeffs,
                    "cov": cov,
                    "m": coeffs[0],
                    "m_err": np.sqrt(cov[0,0]),
                    "b": coeffs[1],
                    "b_err": np.sqrt(cov[1,1]),
                    "R2": R2,
                    "chi2_red": chi2_red,
                    "n_pts": n_pts,
                    "ln_tau_pred": ln_tau_pred,
                    "residuals": residuals,
                    "weights": weights,
                    "inv_T_fine": inv_T_fine,
                    "ln_tau_fine": np.polyval(coeffs,inv_T_fine)
                }
            
            results[e] = {
                "inv_T": inv_T,
                "ln_tau": ln_tau,
                "counts": counts_v,
                "fits": fit_results,
            }

        return results

    @staticmethod
    def calc_activation(tau_data:dict, mode:str="Fit",
                        min_total:int=0) -> dict:
        """
        Calculates the apparent activation energy E_a for each energy
        by fitting ln(tau) vs. 1/T linearly (Arrhenius).

        Parameters
        ----------
        tau_data : dict
            Output of LifetimesIO.gather_tau_vs_T().
        mode : str
            "Fit", "Dir" or "Bin".

        Returns
        -------
        dict mapping e_IS -> {"m": ..., "m_err": ..., "b": ..., "b_err": ...}
        """
        ea_data = {}

        for e, data in tau_data.items():
            sorted_idx = np.argsort(data["temps"])
            temps      = np.array(data["temps"])[sorted_idx]
            inv_T      = 1.0 / temps

            if mode == "Fit":
                tau = np.array(data["tau_ln"])[sorted_idx]
                err = np.array(data["tau_ln_err"])[sorted_idx]
            elif mode == "Dir":
                tau = np.array(data["tau_dir"])[sorted_idx]
                err = np.array(data["tau_dir_err"])[sorted_idx]
            elif mode == "Bin":
                tau = np.array(data["tau_bin"])[sorted_idx]
                err = np.array(data["tau_bin_err"])[sorted_idx]
            else:
                raise ValueError(f"Unknown mode '{mode}'.")
            
            counts_raw = np.array(
                data.get("counts",[1]*len(data["temps"]))
            )[sorted_idx].astype(float)
            total_counts = int(np.sum(counts_raw))

            if min_total > 0 and total_counts < min_total:
                    print(f"INFO: Skipping e={e:.6f} ({mode}): "
                          f"total counts {total_counts} < {min_total}")
                    ea_data[e] = {"m": np.nan, "m_err":np.nan,
                                  "b": np.nan, "b_err":np.nan}
                    continue

            valid = np.isfinite(tau) & (tau > 0) & (temps > 0)

            if np.sum(valid) < 2:
                print(f"WARNING: Not enough valid points for e={e:.4f} ({mode})")
                ea_data[e] = {"m": np.nan, "m_err": np.nan,
                              "b": np.nan, "b_err": np.nan,
                              "R2": np.nan, "chi2": np.nan,
                              "n_points": 0, 
                              "total_counts": total_counts}
                continue

            ln_tau  = np.log(tau[valid])
            weights = np.sqrt(counts_raw[valid])
            n_pts = int(np.sum(valid))

            try:
                coeffs, cov = np.polyfit(inv_T[valid], ln_tau, deg=1,
                                         w=weights,
                                         cov="unscaled")
                
                ln_tau_pred = np.polyval(coeffs,inv_T[valid])
                residuals = ln_tau - ln_tau_pred

                ss_res = np.sum(weights**2 * residuals**2)
                ss_tot = np.sum(weights**2 * (ln_tau - np.average(ln_tau,
                                                                  weights=weights**2))**2)
                
                r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

                dof = n_pts-2
                chi2 = ss_res / dof if dof > 0 else np.nan

                ea_data[e] = {"m":     coeffs[0],
                              "m_err": np.sqrt(cov[0, 0]),
                              "b":     coeffs[1],
                              "b_err": np.sqrt(cov[1, 1]),
                              "R2": r2, "chi2": chi2,
                              "n_points": n_pts, 
                              "total_counts": int(np.sum(counts_raw))}
            except Exception as ex:
                print(f"WARNING: Fit failed for e={e:.4f} ({mode}): {ex}")
                ea_data[e] = {"m": np.nan, "m_err": np.nan,
                              "b": np.nan, "b_err": np.nan,
                              "R2": np.nan, "chi2": np.nan,
                              "n_points": 0, 
                              "total_counts": total_counts}

        return ea_data
    
    @staticmethod
    def calc_quality(path:str,
                     min_counts:list=[1,5,10,20,50],
                     min_totals:list=[0, 50, 100, 200, 500],
                     energies:list=None, thresh:float=1.0, mode:str="Bin") -> dict:
        """
        
        """
        results = {}

        for mc in min_counts:
            tau_data = LifetimesIO.gather_tau_vs_T(
                path=path, energies=energies,thresh=thresh,min_count=mc
            )
            if not tau_data:
                continue
            for mt in min_totals:
                ea = LifetimesFitter.calc_activation(
                    tau_data=tau_data,mode=mode,min_total=mt
                )
                valid_e = [e for e in sorted(ea.keys())
                           if np.isfinite(ea[e]["m"]) and np.isfinite(ea[e]["R2"])]
                
                if not valid_e:
                    results[(mc,mt)] = {
                        "n_valid": 0, 
                        "mean_R2": np.nan, "median_R2": np.nan,
                        "mean_rel_error": np.nan,
                        "mean_chi2": np.nan, "energies": np.array([]),
                        "ea": np.array([]), "ea_err": np.array([]),
                        "R2": np.array([]),
                    }
                    continue

                R2s = np.array([ea[e]["R2"] for e in valid_e])
                ms = np.array([ea[e]["m"] for e in valid_e])
                m_errs = np.array([ea[e]["m_err"] for e in valid_e])
                chi2s = np.array([ea[e]["chi2"] for e in valid_e])
                rel_errs = np.abs(m_errs/ms)

                results[(mc,mt)] = {
                    "n_valid": len(valid_e), 
                    "mean_R2": np.mean(R2s), "median_R2": np.median(R2s),
                    "mean_rel_error": np.mean(rel_errs),
                    "mean_chi2": np.nanmean(chi2s), "energies": np.array(valid_e),
                    "ea": ms, "ea_err": m_errs,
                    "R2": R2s,
                }
        
        return results


# ── LifetimesIO ───────────────────────────────────────────────────────────────

class LifetimesIO:
    """
    Handles all reading and writing for Lifetimes and LifetimesFitter objects.
    """
    @staticmethod
    def read_file(path:str) -> Lifetimes:
        """
        Reads a Lifetimes result file written by write_results().

        Parameters
        ----------
        path : str
        """
        data = np.loadtxt(path)
        life = Lifetimes(
            eis      = data[:, 0].tolist(),
            tau_mean = data[:, 1].tolist(),
            tau_std  = data[:, 2].tolist(),
            dist     = data[:, 3].astype(int).tolist(),
        )
        life.calc_mean()
        return life

    @staticmethod
    def read_tau_vs_T(path:str) -> dict:
        """
        Reads a tau vs. temperature file written by write_tau_vs_T().

        Parameters
        ----------
        path : str
        """
        from collections import defaultdict

        tau_data = defaultdict(lambda: {
            "temps":       [],
            "tau_ln":      [], "tau_ln_err":  [],
            "tau_dir":     [], "tau_dir_err": [],
            "tau_bin":     [], "tau_bin_err": [],
            "counts":      [],
        })

        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                e, temp = float(parts[0]), float(parts[1])
                tau_data[e]["temps"].append(temp)
                tau_data[e]["tau_ln"].append(float(parts[2]))
                tau_data[e]["tau_ln_err"].append(float(parts[3]))
                tau_data[e]["tau_dir"].append(float(parts[4]))
                tau_data[e]["tau_dir_err"].append(float(parts[5]))
                tau_data[e]["tau_bin"].append(float(parts[6]))
                tau_data[e]["tau_bin_err"].append(float(parts[7]))
                if len(parts) > 8:
                    tau_data[e]["counts"].append(int(parts[8]))
                else:
                    tau_data[e]["counts"].append(1)

        return dict(tau_data)

    @staticmethod
    def write_results(life:Lifetimes, path:str) -> None:
        """
        Writes the MB lifetime statistics to a file.

        Parameters
        ----------
        life : Lifetimes
        path : str
        """
        with open(path, "w") as f:
            f.write(f"# <tau>_global = {life.m_tau:.8e}\n")
            for e, tau_m, tau_s, d in zip(
                life.eis, life.tau_mean, life.tau_std, life.dist
            ):
                f.write(f"{e:.10f} {tau_m:.8e} {tau_s:.8e} {d:d}\n")

    @staticmethod
    def write_fits(fitter:LifetimesFitter, path:str, energies:list=None,
                   min_count:int=1) -> None:
        """
        Writes both fit results evaluated at the given energies.

        Parameters
        ----------
        fitter : LifetimesFitter
        path : str
        energies : list, optional
        """
        e_eval, tau_ln, tau_ln_err, tau_dir, tau_dir_err = fitter.calc_tau(
            energies=energies, min_count=min_count
        )

        def _val(arr, i):
            return arr[i] if arr is not None else float("nan")

        with open(path, "w") as f:
            if fitter.fit_ln is not None:
                f.write(f"# ln  fit: a={fitter.fit_ln.a:.10e} "
                        f"b={fitter.fit_ln.b:.10e} c={fitter.fit_ln.c:.10e}\n")
            if fitter.fit_dir is not None:
                f.write(f"# dir fit: a={fitter.fit_dir.a:.10e} "
                        f"b={fitter.fit_dir.b:.10e} c={fitter.fit_dir.c:.10e}\n")
            f.write("# e_IS tau_ln tau_ln_err tau_dir tau_dir_err\n")
            for i, e in enumerate(e_eval):
                f.write(f"{e:.10f} "
                        f"{_val(tau_ln,i):.10e} {_val(tau_ln_err,i):.10e} "
                        f"{_val(tau_dir,i):.10e} {_val(tau_dir_err,i):.10e}\n")

    @staticmethod
    def write_binned(fitter:LifetimesFitter, path:str, energies:list=None,
                     min_count:int=1) -> None:
        """
        Writes the binned tau statistics to a file.

        Parameters
        ----------
        fitter : LifetimesFitter
        path : str
        energies : list, optional
        """
        e_bin, tau_mean, tau_std, dist = fitter.binned_tau(energies=energies,
                                                           min_count=min_count)
        with open(path, "w") as f:
            f.write("# e_IS tau_mean tau_std counts\n")
            for i, e in enumerate(e_bin):
                f.write(f"{e:.10f} {tau_mean[i]:.10e} "
                        f"{tau_std[i]:.10e} {dist[i]:.10e}\n")

    @staticmethod
    def write_tau_vs_T(tau_data:dict, path:str) -> None:
        """
        Writes the tau vs. temperature data to a file.
        Always writes all three variants (ln, dir, bin).

        Parameters
        ----------
        tau_data : dict
            Output of gather_tau_vs_T().
        path : str
        """
        with open(path, "w") as f:
            f.write("# <tau> vs T\n")
            f.write("# e_IS T tau_ln tau_ln_err tau_dir tau_dir_err tau_bin tau_bin_err counts\n")

            for e in sorted(tau_data.keys()):
                data = tau_data[e]
                idx  = np.argsort(data["temps"])

                def _val(key, i):
                    return data[key][i] if key in data else float("nan")

                for i in idx:
                    c = int(data["counts"][i]) if "counts" in data else 0
                    f.write(f"{e:.10f} {data['temps'][i]:.10f} "
                            f"{_val('tau_ln',i):.10e} {_val('tau_ln_err',i):.10e} "
                            f"{_val('tau_dir',i):.10e} {_val('tau_dir_err',i):.10e} "
                            f"{_val('tau_bin',i):.10e} {_val('tau_bin_err',i):.10e} "
                            f"{c:d}\n")

    @staticmethod
    def gather_tau_vs_T(path:str, energies:list=None, thresh:float=1.0,
                        min_count:int=1) -> None:
        """
        Calculates mean tau at different temperatures for given energies,
        always computing all three variants (ln fit, dir fit, binning).

        Parameters
        ----------
        path : str
            Glob pattern matching Lifetimes result files.
        energies : list, optional
        thresh : float
        min_count : int
            Minimum occupation for a bin to be considered.
        """
        import re
        from collections import defaultdict

        def extract_temp(filepath):
            match = re.search(r'T(\d+\.\d+)', filepath)
            if match:
                return float(match.group(1))
            raise ValueError(f"No temperature found in: {filepath}")

        files = glob.glob(path)

        if energies is None:
            e_min, e_max = np.inf, -np.inf
            for file in files:
                try:
                    life_tmp = LifetimesIO.read_file(file)
                    if life_tmp.filter_tau(thresh=thresh):
                        e_min = min(e_min,min(life_tmp.eis))
                        e_max = max(e_max,max(life_tmp.eis))
                except Exception as e:
                    print(f"WARNING: Skipping {file}: {e}")
            if np.isinf(e_min): 
                print("ERROR: No valid files found")
                return {}
            bin_edges = np.linspace(e_min,e_max,51)
        else:
            bin_edges = np.asarray(energies)

        if len(bin_edges) < 2:
            print("ERROR: Need at least 2 edge values.")
            return {}
        
        bin_centers = ((bin_edges[:-1] + bin_edges[1:]) / 2).tolist()

        tau_data = defaultdict(lambda: {
            "temps":       [],
            "tau_ln":      [], "tau_ln_err":  [],
            "tau_dir":     [], "tau_dir_err": [],
            "tau_bin":     [], "tau_bin_err": [],
            "counts":      [],
        })

        for file in files:
            temp = extract_temp(file)
            life = LifetimesIO.read_file(file)

            if not life.filter_tau(thresh=thresh):
                print(f"WARNING: Skipping T={temp:.4f}, filtering failed.")
                continue

            fitter = LifetimesFitter(life)
            fitter.binned_tau(energies=bin_edges.tolist(),min_count=min_count)
            fitter.fit_ln_tau()
            fitter.fit_dir_tau(min_count=min_count)
            e_eval, tau_ln, tau_ln_err, tau_dir, tau_dir_err = fitter.calc_tau(
                min_count=min_count
            )
            
            tau_bin = fitter.binned_tau_mean
            tau_bin_err = fitter.binned_tau_std
            bin_counts = fitter.binned_dist

            for i, e in enumerate(e_eval):
                tau_data[e]["temps"].append(temp)
                tau_data[e]["tau_ln"].append(
                    tau_ln[i] if tau_ln is not None else np.nan)
                tau_data[e]["tau_ln_err"].append(
                    tau_ln_err[i] if tau_ln_err is not None else np.nan)
                tau_data[e]["tau_dir"].append(
                    tau_dir[i] if tau_dir is not None else np.nan)
                tau_data[e]["tau_dir_err"].append(
                    tau_dir_err[i] if tau_dir_err is not None else np.nan)
                tau_data[e]["tau_bin"].append(
                    tau_bin[i] if tau_bin is not None else np.nan)
                tau_data[e]["tau_bin_err"].append(
                    tau_bin_err[i] if tau_bin_err is not None else np.nan)
                tau_data[e]["counts"].append(
                    int(bin_counts[i]) if bin_counts is not None else 0)

        return dict(tau_data)
    
    @staticmethod
    def check_tauT(path:str, energies:list=None, thresh:float=0.1,
                   min_count:int=1,
                   mask_ln:bool=True) -> dict:
        """
        
        """
        import re
        from collections import defaultdict

        def extract_temp(filepath):
            match = re.search(r'T(\d+\.\d+)',filepath)
            if match:
                return float(match.group(1))
            raise ValueError(f"No temperature found in: {filepath}")
        
        files = glob.glob(path)

        if energies is None:
            e_min, e_max = np.inf, -np.inf
            for file in files:
                try:
                    life_tmp = LifetimesIO.read_file(file)
                    e_min = min(e_min,min(life_tmp.eis))
                    e_max = max(e_max,max(life_tmp.eis))
                except Exception as ex:
                    print(f"WARNING: Skipping {file}: {ex}")
            bin_edges = np.linspace(e_min,e_max,51)
        else:
            bin_edges = np.asarray(energies)

        bin_centers = ((bin_edges[:-1] + bin_edges[1:]) / 2).tolist()
        n_bins = len(bin_centers)

        tau_data = defaultdict(lambda: {
            "temps": [],
            "tau_ln": [], "tau_ln_err": [],
            "tau_dir": [], "tau_dir_err": [],
            "tau_bin": [], "tau_bin_err": [],
            "counts": [],
        })

        for file in files:
            try:
                temp = extract_temp(file)
            except ValueError as ex:
                print(f"WARNGING: {ex}")
                continue

            life = LifetimesIO.read_file(file)
            #LifetimesPlotter.check_plot_binned(life,nr_bins_list=[0],
            #                                   min_count=min_count,energies=bin_edges.tolist(),
            #                                   save=str(temp)+"_binned")
            #LifetimesPlotter.check_plot_fitted(life,nr_bins_list=[0],
            #                                   min_count=min_count,energies=bin_edges.tolist(),
            #                                   save=str(temp)+"_fitted")
            fitter = LifetimesFitter(life)

            diag = fitter.check_fitted_tau(min_count=min_count,
                                           energies=bin_edges.tolist())
            stats = diag["bin_stats"]
            fits = diag["fits"]

            tau_bin = stats["tau_mean"]
            tau_bin_sem = stats["tau_sem"]
            counts = stats["counts"]

            fit_ln = fits.get("inv_sig")
            if fit_ln is not None:
                centers_arr = stats["centers"]
                coeffs = fit_ln["coeffs"]
                cov = fit_ln["cov"]
                ln_tau_eval = np.polyval(coeffs,centers_arr)
                tau_ln_eval = np.exp(ln_tau_eval)
                tau_ln_err_eval = np.array([
                    tau_ln_eval[i] * np.sqrt(
                        np.array([centers_arr[i]**2,centers_arr[i], 1])
                        @ cov
                        @ np.array([centers_arr[i]**2,centers_arr[i],1])
                    )
                    for i in range(n_bins)
                ])
            else:
                tau_ln_eval = np.full(n_bins,np.nan)
                tau_ln_err_eval = np.full(n_bins,np.nan)

            if mask_ln:
                tau_ln_eval[counts < min_count] = np.nan
                tau_ln_err_eval[counts < min_count] = np.nan

            for i,e in enumerate(bin_centers):
                tau_data[e]["temps"].append(temp)
                tau_data[e]["tau_ln"].append(float(tau_ln_eval[i]))
                tau_data[e]["tau_ln_err"].append(float(tau_ln_err_eval[i]))
                tau_data[e]["tau_bin"].append(float(tau_bin[i]))
                tau_data[e]["tau_bin_err"].append(float(tau_bin_sem[i]))
                tau_data[e]["counts"].append(float(counts[i]))

        return dict(tau_data)

    @staticmethod
    def write_counts(tau_data:dict,path:str,energies:list=None,
                     min_total:int=0) -> None:
        """
        
        """
        if energies is not None:
            available = sorted(tau_data.keys())
            energies_sorted = []
            for e in energies:
                nearest = min(available,key=lambda x: abs(x-e))
                if nearest not in energies_sorted:
                    energies_sorted.append(nearest)
            energies_sorted = sorted(energies_sorted)
        else:
            energies_sorted = sorted(tau_data.keys())

        with open(path,"w") as f:
            f.write("# e_IS total_counts\n")
            for e in energies_sorted:
                counts = tau_data[e].get("counts",[])
                total = sum(int(c) for c in counts if np.isfinite(c))
                if total < min_total:
                    continue
                f.write(f"{e:.10f} {total:d}\n")
        
        return

    @staticmethod
    def write_activ(activ_data:dict,path:str,mode:str="inv_sig",energies:list=None):
        """
        
        """
        available = sorted(activ_data.keys())
        if energies is not None:
            energies_out = sorted(set(
                min(available, key=lambda x: abs(x - e)) for e in energies
            ))
        else:
            energies_out = available

        with open(path, "w") as f:
            f.write(f"# E_a vs e_IS  (mode={mode})\n")
            f.write("# e_IS m m_err b b_err R2 chi2_red n_pts\n")
            for e in energies_out:
                entry = activ_data[e]
                if entry is None:
                    continue
                fit = entry["fits"].get(mode)
                if fit is None:
                    continue
                f.write(f"{e:.10f} {fit['m']:.8e} {fit['m_err']:.8e} "
                        f"{fit['b']:.8e} {fit['b_err']:.8e} "
                        f"{fit['R2']:.6f} {fit['chi2_red']:.6f} "
                        f"{fit['n_pts']:d}\n")
        return


# ── LifetimesPlotter ──────────────────────────────────────────────────────────

class LifetimesPlotter:
    """
    Handles all plotting for Lifetimes and LifetimesFitter objects.
    Not intended for use on HPC clusters – imports matplotlib/seaborn/scienceplots.
    """
    def __init__(self, life:Lifetimes, fitter:LifetimesFitter=None):
        self.life   = life
        self.fitter = fitter
        self.plt    = self._setup_plotting()

    @staticmethod
    def _setup_plotting():
        import seaborn as sns
        import matplotlib as mpl
        mpl.use('TkAgg')
        import scienceplots
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'sans-serif'
        plt.style.use(['science', 'nature'])
        return plt

    @staticmethod
    def _save(fig, name):
        fig.savefig(f"{name}.png")
        fig.savefig(f"{name}.pdf")

    def plot_dist(self, save:str="dist") -> None:
        """
        Plots tau(e_IS) and the MB distribution phi(e_IS).
        Overlays fit curves and binned data if available via self.fitter.
        """
        fig, (ax1, ax2) = self.plt.subplots(1, 2, figsize=(3.39, 2.0))

        ax1.plot(self.life.eis, self.life.tau_mean,
                 marker="o", markersize=4, linestyle="None")
        ax1.axhline(self.life.m_tau, linestyle="--")

        if self.fitter is not None:
            e_fit = np.linspace(min(self.life.eis), max(self.life.eis), 200)

            if self.fitter.fit_ln is not None:
                tau_fit, _ = self.fitter.fit_ln.eval(e_fit)
                ax1.plot(e_fit, tau_fit, "--", label=r"$\ln\tau$ fit")

            if self.fitter.fit_dir is not None:
                tau_fit, _ = self.fitter.fit_dir.eval(e_fit)
                ax1.plot(e_fit, tau_fit, "--", label=r"$\tau$ fit")

            if self.fitter.binned_energies is not None:
                ax1.plot(self.fitter.binned_energies, self.fitter.binned_tau_mean,
                         marker="x", markersize=5, linestyle="None", label="binned")

        ax1.set_yscale("log")
        ax1.set_ylim(min(self.life.tau_mean), max(self.life.tau_mean))
        ax1.set_xlabel(r"$e_{IS}$")
        ax1.set_ylabel(r"$\langle \tau (e_{IS}) \rangle$")
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.legend()

        ax2.plot(self.life.eis, self.life.dist,
                 marker="o", markersize=4, linestyle="None")
        ax2.set_xlabel(r"$e_{IS}$")
        ax2.set_ylabel(r"$\varphi (e_{IS})$")
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')

        self.plt.tight_layout()
        self._save(fig, save)
        self.plt.show(block=False)
        plt.pause(3)
        self.plt.close()
    
    @staticmethod
    def _filter_by_min_total(tau_data:dict, energies:list, min_total:int) -> list:
        """
        
        """
        if min_total <= 0:
            return energies
        filtered = []
        for e in energies:
            counts = tau_data[e].get("counts", [])
            total = sum(c for c in counts if np.isfinite(c))
            if total >= min_total:
                filtered.append(e)
            else:
                print(f"INFO: Omitting e={e:.6f} from plot "
                      f"(total counts {total:.0f} < {min_total})")
        return np.array(filtered)

    @staticmethod
    def plot_tau_vs_T(tau_data:dict, save:str="tau_vs_T",energies:list=None,
                      min_total:int=0, 
                      ea_ln:dict=None, ea_dir:dict=None, ea_bin:dict=None) -> None:
        """
        Plots tau(e_IS, T) vs. 1/T for all energies in tau_data.
        Always plots all three variants (ln, dir, bin) as separate panels.

        Parameters
        ----------
        tau_data : dict
            Output of LifetimesIO.gather_tau_vs_T().
        save : str
            Base name for output files.
        """
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        plt = LifetimesPlotter._setup_plotting()

        if energies is not None:
            available = sorted(tau_data.keys())
            energies_sorted = []
            for e in energies:
                nearest = min(available,key=lambda x: abs(x-e))
                if nearest not in energies_sorted:
                    energies_sorted.append(nearest)
            energies_sorted = sorted(energies_sorted)
        else:
            energies_sorted = sorted(tau_data.keys())

        energies_sorted = LifetimesPlotter._filter_by_min_total(
            tau_data=tau_data,energies=energies_sorted,min_total=min_total
        )

        colors = plt.cm.viridis(np.linspace(0, 1, len(energies_sorted)))

        panels = [
            ("tau_ln",  "tau_ln_err",  r"$\ln\tau$ fit", ea_ln),
            ("tau_dir", "tau_dir_err", r"$\tau$ fit", ea_dir),
            ("tau_bin", "tau_bin_err", "binned", ea_bin),
        ]

        fig, axes = plt.subplots(1, 3, figsize=(3*3.39, 2.0))

        for ax, (tau_key, err_key, title, ea_data) in zip(axes, panels):
            all_inv_T = []
            for i, e in enumerate(energies_sorted):
                data  = tau_data[e]
                idx   = np.argsort(data["temps"])
                temps = np.array(data["temps"])[idx]
                taus  = np.array(data[tau_key])[idx]
                errs  = np.array(data[err_key])[idx]
                valid = ~np.isnan(taus)

                inv_T = 1.0 / temps[valid]
                all_inv_T.extend(inv_T.tolist())

                ax.errorbar(1/temps[valid], taus[valid], yerr=errs[valid],
                            marker="o", markersize=4, linestyle="--",
                            capsize=3, capthick=1, color=colors[i])
                
                if ea_data is not None and e in ea_data:
                    m = ea_data[e]["m"]
                    b = ea_data[e]["b"]
                    if np.isfinite(m) and np.isfinite(b) and len(inv_T) >= 2:
                        inv_T_fit = np.linspace(inv_T.min(), inv_T.max(), 50)
                        tau_fit = np.exp(m*inv_T_fit+b)
                        ax.plot(inv_T_fit,tau_fit,linestyle="--",
                                linewidth=0.7,color=colors[i],alpha=0.7)

            ax.set_yscale("log")
            ax.set_title(title)
            ax.set_xlabel(r"$1/T$")
            ax.set_ylabel(r"$\langle \tau (e_{IS},T) \rangle$")
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

        norm = Normalize(vmin=min(energies_sorted), vmax=max(energies_sorted))
        sm   = ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        divider = make_axes_locatable(axes[-1])
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        cbar    = fig.colorbar(sm, cax=cax)
        cbar.set_label(r"$e_{IS}$")

        plt.tight_layout()
        LifetimesPlotter._save(fig, save)
        plt.show(block=False)
        plt.pause(3)
        plt.close()

    @staticmethod
    def plot_activation(ea_ln:dict, ea_dir:dict, ea_bin:dict, 
                        save:str="activation",energies:list=None,
                        tau_data:dict=None,min_total:int=0) -> None:
        """
        Plots E_a(e_IS) for all three variants (ln, dir, bin) as separate panels.

        Parameters
        ----------
        ea_ln  : dict
            Output of LifetimesFitter.calc_activation(tau_data, mode="Fit").
        ea_dir : dict
            Output of LifetimesFitter.calc_activation(tau_data, mode="Dir").
        ea_bin : dict
            Output of LifetimesFitter.calc_activation(tau_data, mode="Bin").
        save : str
            Base name for output files.
        """
        plt = LifetimesPlotter._setup_plotting()
        fig, axes = plt.subplots(1, 3, figsize=(3*3.39, 2.0))

        for ax, ea_data, title in [
            (axes[0], ea_ln,  r"$\ln\tau$ fit"),
            (axes[1], ea_dir, r"$\tau$ fit"),
            (axes[2], ea_bin, "binned"),
        ]:
            if energies is not None:
                available = sorted(ea_data.keys())
                energies_plot = np.array(sorted(
                    set(min(available,key=lambda x: abs(x-e)) for e in energies)
                ))
            else:
                energies_plot = np.array(sorted(ea_data.keys()))

            if tau_data is not None:
                energies_plot = LifetimesPlotter._filter_by_min_total(
                    tau_data=tau_data,energies=energies_plot,min_total=min_total
                )

            ea       = np.array([ea_data[e]["m"]     for e in energies_plot])
            ea_err   = np.array([ea_data[e]["m_err"] for e in energies_plot])
            tau0     = np.array([ea_data[e]["b"]     for e in energies_plot])
            tau0_err = np.array([ea_data[e]["b_err"]     for e in energies_plot])
            valid    = np.isfinite(ea)

            ax.errorbar(energies_plot[valid], ea[valid], yerr=ea_err[valid],
                        marker="o", markersize=4, linestyle="None",
                        capsize=3, capthick=1)

            ax.set_title(title)
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$E_a (e_{IS})$")
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

            #inset = ax.inset_axes([0.55,0.55,0.42,0.42])
            #inset.errorbar(energies_plot[valid],tau0[valid],yerr=tau0_err[valid],
            #               marker="o",markersize=4,linestyle="None",
            #               capsize=3,capthick=1)
            #inset.set_xlabel(r"$e_{IS}$",fontsize=6)
            #inset.set_ylabel(r"$\ln\tau_0$",fontsize=6)
            #inset.xaxis.set_ticks_position('bottom')
            #inset.yaxis.set_ticks_position('left')
            #inset.tick_params(labelsize=5)

        plt.tight_layout()
        LifetimesPlotter._save(fig, save)
        plt.show(block=False)
        plt.pause(3)
        plt.close()

    @staticmethod
    def plot_global_tau(path:str, save:str="global_tau") -> None:
        """
        Plots the global mean tau vs. T and ln(tau) vs. 1/T
        from a set of Lifetimes result files.

        Parameters
        ----------
        path : str
            Glob pattern matching Lifetimes result files.
        save : str
            Base name for output files.
        """
        import re

        plt = LifetimesPlotter._setup_plotting()

        def extract_temperature(filepath):
            match = re.search(r'T(\d+\.\d+)', filepath)
            if match:
                return float(match.group(1))
            raise ValueError(f"No temperature found in: {filepath}")

        files = sorted(glob.glob(path))
        temps, tau_means, tau_stds = [], [], []

        for file in files:
            try:
                life = LifetimesIO.read_file(file)
                temps.append(extract_temperature(file))
                tau_means.append(life.m_tau)
                tau_stds.append(np.std(life.tau_mean))
            except Exception as e:
                print(f"WARNING: Skipping {file}: {e}")

        temps     = np.array(temps)
        tau_means = np.array(tau_means)
        tau_stds  = np.array(tau_stds)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(3.39, 2.0))

        ax1.errorbar(temps, tau_means, yerr=tau_stds,
                     marker="o", markersize=4, linestyle="None",
                     capsize=3, capthick=1)
        ax1.set_xlabel(r"$T$")
        ax1.set_ylabel(r"$\langle \tau \rangle$")
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')

        tau_std_log = tau_stds / tau_means
        ax2.errorbar(1/temps, np.log(tau_means), yerr=tau_std_log,
                     marker="o", markersize=4, linestyle="None",
                     capsize=3, capthick=1)
        ax2.set_xlabel(r"$1/T$")
        ax2.set_ylabel(r"$\ln \langle \tau \rangle$")
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')

        plt.tight_layout()
        LifetimesPlotter._save(fig, save)
        plt.show(block=False)
        plt.pause(3)
        plt.close()

    @staticmethod
    def plot_debug(life: Lifetimes, fitter: LifetimesFitter,
                   save: str = "debug") -> None:
        """
        Debug plots for fit_ln_tau, fit_dir_tau, and binned_tau.

        Panel 1 – ln tau fit:  log plot + linear inset
        Panel 2 – dir tau fit: log plot + linear inset + bin means
        Panel 3 – binned:      bin means with edges + counts at top
        """
        plt = LifetimesPlotter._setup_plotting()
        fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.5))

        e_raw  = np.array(life.eis)
        t_raw  = np.array(life.tau_mean)
        e_fine = np.linspace(e_raw.min(), e_raw.max(), 300)

        # ── Panel 1: ln tau fit ───────────────────────────────────────────────
        ax = axes[0]
        ax.plot(e_raw, t_raw, marker="o", markersize=3,
                linestyle="None",alpha=0.5,
                markeredgecolor="none")

        if fitter.fit_ln is not None:
            tau_fit, tau_err = fitter.fit_ln.eval(e_fine)
            ax.plot(e_fine, tau_fit, lw=1.2, color="red",linestyle="--")

        ax.set_yscale("log")
        ax.set_xlabel(r"$e_{IS}$")
        ax.set_ylabel(r"$\langle\tau\rangle$")
        ax.set_title(r"$\ln\tau$ fit")
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")

        ins = ax.inset_axes([0.55,0.55,0.42,0.38])
        ins.plot(e_raw,t_raw,marker="o",markersize=3,
                 linestyle="None",alpha=0.5,
                 markeredgecolor="none")
        ins.plot(e_fine,tau_fit,lw=1.0,color="red")
        ins.set_xlabel(r"$e_{IS}$",fontsize=5)
        ins.set_ylabel(r"$\tau$",fontsize=5)
        ins.tick_params(labelsize=4)
        ins.xaxis.set_ticks_position("bottom")
        ins.yaxis.set_ticks_position("left")

        # ── Panel 2: dir tau fit ──────────────────────────────────────────────
        ax = axes[1]
        ax.plot(e_raw, t_raw, marker="o", markersize=3,
                linestyle="None",alpha=0.5,
                markeredgecolor="none")

        if fitter.fit_dir is not None:
            tau_fit, tau_err = fitter.fit_dir.eval(e_fine)
            ax.plot(e_fine, tau_fit, lw=1.2,color="red",linestyle="--")

        ax.set_yscale("log")
        ax.set_xlabel(r"$e_{IS}$")
        ax.set_ylabel(r"$\langle\tau\rangle$")
        ax.set_title(r"$\tau$ fit")
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")

        if fitter.binned_energies is not None:
            valid = np.isfinite(fitter.binned_tau_mean)
            ax.plot(fitter.binned_energies[valid],
                    fitter.binned_tau_mean[valid],
                    marker="x", markersize=5, linestyle="None",
                    color="red", label="bin means", zorder=4)

            ins = ax.inset_axes([0.55,0.55,0.42,0.38])
            ins.plot(e_raw,t_raw,marker="o",markersize=3,
                     linestyle="None",alpha=0.5,
                     markeredgecolor="none")
            ins.plot(e_fine,tau_fit,lw=1.0,color="red")
            ins.set_xlabel(r"$e_{IS}$",fontsize=5)
            ins.set_ylabel(r"$\tau$",fontsize=5)
            ins.tick_params(labelsize=4)
            ins.xaxis.set_ticks_position("bottom")
            ins.yaxis.set_ticks_position("left")

        # ── Panel 3: binned ───────────────────────────────────────────────────────
        ax = axes[2]

        if fitter.binned_energies is not None:
            be    = fitter.binned_energies
            bm    = fitter.binned_tau_mean
            bs    = fitter.binned_tau_std
            bd    = fitter.binned_dist
            valid = np.isfinite(bm) & (bd > 0)

            # Bin-Edges rekonstruieren
            half  = (be[1] - be[0]) / 2 if len(be) > 1 else 0
            edges = np.concatenate([[be[0] - half],
                                     (be[:-1] + be[1:]) / 2,
                                     [be[-1] + half]])

            # Bin-Mittelwerte mit Fehlerbalken
            ax.errorbar(be[valid], bm[valid],
                        yerr=bs[valid],
                        marker="o", markersize=3, linestyle="None",
                        capsize=2, capthick=0.8, color="steelblue")

            # Bin-Grenzen als vertikale Linien
            for edge in edges:
                ax.axvline(edge, color="gray", lw=0.5, ls="--", alpha=0.5)

            ax.set_yscale("log")
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$\langle\tau\rangle$")
            ax.set_title("binned")
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")

            # Inset: Histogramm der Bin-Besetzung
            ins = ax.inset_axes([0.55, 0.55, 0.42, 0.38])
            ins.bar(be, bd, width=2*half, color="steelblue", alpha=0.7,
                    edgecolor="white", linewidth=0.3)
            ins.set_xlabel(r"$e_{IS}$", fontsize=5)
            ins.set_ylabel(r"$N$",      fontsize=5)
            ins.tick_params(labelsize=4)
            ins.xaxis.set_ticks_position("bottom")
            ins.yaxis.set_ticks_position("left")

        plt.tight_layout()
        LifetimesPlotter._save(fig, save)
        plt.show(block=False)
        plt.pause(3)
        plt.close()

    @staticmethod
    def plot_counts(tau_data:dict, save:str = "counts", energies:list=None,
                    min_total:int=0, normalize:str="energy") -> None:
        """
        
        """
        from matplotlib.colors import LogNorm
        plt = LifetimesPlotter._setup_plotting()

        if energies is not None:
            available = sorted(tau_data.keys())
            energies_sorted = []
            for e in energies:
                nearest = min(available,key=lambda x: abs(x-e))
                if nearest not in energies_sorted:
                    energies_sorted.append(nearest)
            energies_sorted = sorted(energies_sorted)
        else:
            energies_sorted = sorted(tau_data.keys())

        energies_sorted = LifetimesPlotter._filter_by_min_total(
            tau_data=tau_data, energies=energies_sorted,
            min_total=min_total
        )

        all_temps = sorted(set(
            t for e in energies_sorted for t in tau_data[e]["temps"]
        ))

        dist_matrix = np.zeros((len(energies_sorted),len(all_temps)))

        for i,e in enumerate(energies_sorted):
            data = tau_data[e]
            for j,t in enumerate(all_temps):
                if t in data["temps"]:
                    idx = data["temps"].index(t)
                    dist_matrix[i,j] = data["counts"][idx]
        
        counts_per_energy = dist_matrix.sum(axis=1)
        counts_per_temp = dist_matrix.sum(axis=0)

        if normalize == "temp":
            plot_matrix = np.zeros_like(dist_matrix)
            for j in range(len(all_temps)):
                col_sum = counts_per_temp[j]
                if col_sum > 0:
                    plot_matrix[:,j] = dist_matrix[:,j] / col_sum
            clabel = r"$p(e_{IS} | T)$"
            use_log = False
        elif normalize == "energy":
            plot_matrix = np.zeros_like(dist_matrix)
            for j in range(len(energies_sorted)):
                col_sum = counts_per_energy[j]
                if col_sum > 0:
                    plot_matrix[j,:] = dist_matrix[j,:] / col_sum
            clabel = r"$p(T | e_{IS})$"
            use_log = False
        else:
            plot_matrix = dist_matrix
            clabel = r"$N$"
            use_log = True

        fig = plt.figure(figsize=(3.39*2,3.39*1.5))
        gs = fig.add_gridspec(2,2,width_ratios=[4,1],height_ratios=[1,4],
                              hspace=0.05,wspace=0.05)
        ax_heat = fig.add_subplot(gs[1,0])
        ax_top = fig.add_subplot(gs[0,0],sharex=ax_heat)
        ax_right = fig.add_subplot(gs[1,1],sharey=ax_heat)

        temp_labels = [f"{t:.3f}" for t in all_temps]
        energy_labels = [f"{e:.4f}" for e in energies_sorted]

        if use_log == False:
            im = ax_heat.imshow(plot_matrix,aspect="auto",origin="lower",
                                cmap="viridis",
                                extent=[-0.5,len(all_temps)-0.5,
                                        -0.5,len(energies_sorted)-0.5])
        else:
            vmin = max(1,plot_matrix[plot_matrix > 0].min()) if np.any(plot_matrix > 0) else 1
            im = ax_heat.imshow(plot_matrix,aspect="auto",origin="lower",
                                cmap="viridis",
                                norm=LogNorm(vmin=vmin,vmax=plot_matrix.max()),
                                extent=[-0.5,len(all_temps)-0.5,
                                        -0.5,len(energies_sorted)-0.5])

        ax_heat.set_xticks(range(len(all_temps)))
        ax_heat.set_xticklabels(temp_labels,rotation=45,ha="right",fontsize=5)
        ax_heat.set_yticks(range(len(energies_sorted)))
        ax_heat.set_yticklabels(energy_labels,fontsize=5)
        ax_heat.set_xlabel(r"$T$")
        ax_heat.set_ylabel(r"$e_{IS}$")

        ax_top.bar(range(len(all_temps)),counts_per_temp,
                   alpha=0.7)
        ax_top.set_ylabel(r"$\sum N$",fontsize=6)
        ax_top.tick_params(labelleft=False,labelsize=5)
        ax_top.xaxis.set_ticks_position('bottom')
        ax_top.yaxis.set_ticks_position('left')
        ax_top.set_yscale("log")

        ax_right.barh(range(len(energies_sorted)),counts_per_energy,
                      alpha=0.7)
        ax_right.set_ylabel(r"$\sum N$",fontsize=6)
        ax_right.tick_params(labelleft=False,labelsize=5)
        ax_right.xaxis.set_ticks_position('bottom')
        ax_right.yaxis.set_ticks_position('left')

        cax = fig.add_axes([0.92,0.12,0.02,0.55])
        cbar = fig.colorbar(im,cax=cax)
        cbar.set_label(clabel,fontsize=6)
        cbar.ax.tick_params(labelsize=5)

        plt.tight_layout()
        LifetimesPlotter._save(fig,save)
        plt.show(block=False)
        plt.pause(3)
        plt.close()

    @staticmethod
    def check_plot_binned(life:Lifetimes,nr_bins_list:list=[25,50,100],
                          min_count:int=1,save:str="binning",
                          energies:list=None):
        """
        
        """
        plt = LifetimesPlotter._setup_plotting()

        fitter = LifetimesFitter(life=life)
        e_raw = np.array(life.eis)
        t_raw = np.array(life.tau_mean)
        d_raw = np.array(life.dist)

        all_stats = {}

        for nb in nr_bins_list:
            stats = fitter.check_binned_tau(
                nr_bins=nb,min_count=min_count,energies=energies
            )
            all_stats[nb] = stats

            bins = len(stats["centers"])

            fig,axes = plt.subplots(3,1,figsize=(7,5.5))
            fig.suptitle(
                rf"$n_{{\mathrm{{bins}}}} = {bins}$, "
                rf"$N_{{\mathrm{{min}}}} = {min_count}$",
            )

            valid = np.isfinite(stats["tau_mean"])

            ax = axes[0]
            sizes = np.clip(d_raw*2,2,50)
            ax.scatter(e_raw,t_raw,s=sizes,alpha=0.3,edgecolors="none")
            ax.errorbar(
                stats["centers"][valid],
                stats["tau_mean"][valid],
                yerr = stats["tau_sem"][valid],
                marker = "s", markersize = 4, linestyle = "None",
                color = "red", label = r"$\bar\tau\pm$ SEM"
            )
            for edge in stats["edges"]:
                ax.axvline(edge,color="gray",lw=0.3,ls="--",alpha=0.3)
            ax.set_yscale("log")
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$\tau$")

            ax = axes[1]
            width = stats["edges"][1] - stats["edges"][0]
            ax.bar(stats["centers"],stats["counts"],
                   width = width*0.9)
            ax.axhline(min_count,color="red",ls="--")
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$N$ (counts)")
            ax.set_yscale("log")

            ax = axes[2]
            ax.plot(stats["centers"][valid],
                    stats["rel_sem"][valid],
                    marker="o",markersize=3,linestyle="-")
            ax.axhline(0.5,color="orange",ls="--",lw=0.8,
                       label="rel SEM=0.5")
            ax.axhline(0.1,color="green",ls="--",lw=0.8,
                       label="rel SEM=0.1")
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"SEM / $\bar\tau$")
            ax.set_yscale("log")
            
            plt.tight_layout()
            LifetimesPlotter._save(fig,save+str(bins))
            plt.show(block=False)
            plt.pause(3)
            plt.close()
    
    @staticmethod
    def check_plot_fitted(life,nr_bins_list:list=[25,50,100],
                          min_count:int=1,save:str="fit_check",
                          energies:list=None):
        """
        
        """
        plt = LifetimesPlotter._setup_plotting()
        fitter = LifetimesFitter(life)

        e_raw = np.array(life.eis)
        t_raw = np.array(life.tau_mean)

        for nb in nr_bins_list:
            diag = fitter.check_fitted_tau(nr_bins=nb,min_count=min_count,energies=energies)
            if diag is None:
                continue
            stats = diag["bin_stats"]
            fits = diag["fits"]
            valid_weights = {k:v for k,v in fits.items() if v is not None}
            n_w = len(valid_weights)
            bins = len(stats["centers"])
            if n_w == 0:
                continue

            fig,axes = plt.subplots(3,n_w,figsize=(3.5*n_w,7),squeeze=False)

            fig.suptitle(
                rf"$n_{{\mathrm{{bins}}}} = {bins}$, "
                rf"$N_{{\mathrm{{min}}}} = {min_count}$",
            )

            for col,(wname,f) in enumerate(
                valid_weights.items()
            ):
                label = {
                    "sqrt_n": r"$w = \sqrt{N}$",
                    "inv_sig": r"$w = 1/\sigma(\ln\tau)$",
                    "uniform": r"$w = 1$"
                }.get(wname,wname)

                ax = axes[0,col]
                ax.scatter(e_raw,t_raw,s=1,alpha=0.5,color="gray")
                valid = np.isfinite(stats["tau_mean"])
                ax.errorbar(
                    stats["centers"][valid],
                    stats["tau_mean"][valid],
                    yerr=stats["tau_sem"][valid],
                    marker="s", markersize=3,
                    linestyle="None", color="steelblue",
                    capsize=1, capthick=0.5, zorder=2
                )
                ax.plot(f["e_fine"],f["tau_fine"],color="red",lw=1.2,zorder=3)
                ax.set_yscale("log")
                ax.set_xlabel(r"$e_{IS}$")
                ax.set_ylabel(r"$\tau$")
                ax.set_title(
                    f"{label}\n"
                    rf"$R^2 = {f['R2']:.3f}$, "
                    rf"$\chi^2_{{\mathrm{{red}}}} = {f['chi2_red']:.1f}$"
                )

                ax = axes[1,col]
                ax.axhline(0,color="gray")
                colors = plt.cm.viridis(
                    f["counts"] / f["counts"].max()
                )
                ax.scatter(f["e_valid"],f["residuals"],c=colors,s=15,zorder=2)

                if wname == "inv_sig":
                    ax.fill_between(
                        f["e_valid"],
                        -f["sigma_ln"],
                        +f["sigma_ln"],
                        alpha = 0.2, color = "red"
                    )

                ax.set_xlabel(r"$e_{IS}$")
                ax.set_ylabel(
                    r"$\ln\tau - \ln\tau{\mathrm{fit}}$"
                )
                ax.set_title("Residues")

                ax = axes[2,col]
                ax.scatter(f["counts"],np.abs(f["residuals"]),s=10,alpha=0.6)
                ax.set_xlabel(r"$N$ (counts)")
                ax.set_ylabel(r"$| \mathrm{Residuum} |$")
                ax.set_xscale("log")
            
            plt.tight_layout()
            LifetimesPlotter._save(fig,save+str(bins))
            plt.show(block=False)
            plt.pause(3)
            plt.close()

    @staticmethod
    def check_plot_tauT(tau_data:dict,save:str="check_tauT",
                        min_total:int=1,energies:list=None,
                        ea_ln:dict=None,ea_bin:dict=None):
        """
        
            """
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize, LogNorm
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        plt = LifetimesPlotter._setup_plotting()

        if energies is not None:
            available       = sorted(tau_data.keys())
            energies_sorted = sorted(set(
                min(available, key=lambda x: abs(x - e)) for e in energies
            ))
        else:
            energies_sorted = sorted(tau_data.keys())

        if min_total > 0:
            energies_sorted = [
                e for e in energies_sorted
                if sum(
                    c for c in tau_data[e].get("counts",[])
                    if np.isfinite(c)
                ) >= min_total
            ]

        all_temps = sorted(set(
            t for e in energies_sorted for t in tau_data[e]["temps"]
        ))

        colors = plt.cm.viridis(np.linspace(0, 1, len(energies_sorted)))

        fig, axes = plt.subplots(1, 3, figsize=(3 * 3.39, 2.0))

        # ── Panel 1: tau_ln vs 1/T ────────────────────────────────────────────────
        ax = axes[0]
        for i, e in enumerate(energies_sorted):
            data  = tau_data[e]
            idx   = np.argsort(data["temps"])
            temps = np.array(data["temps"])[idx]
            taus  = np.array(data["tau_ln"])[idx]
            errs  = np.array(data["tau_ln_err"])[idx]

            valid = np.isfinite(taus) & (taus > 0)
            if not np.any(valid):
                continue

            inv_T = 1.0 / temps[valid]
            ax.errorbar(inv_T, taus[valid], yerr=errs[valid],
                        marker="o", markersize=4, linestyle="--",
                        capsize=3, capthick=1, color=colors[i])

            if ea_ln is not None and e in ea_ln:
                m = ea_ln[e]["m"]
                b = ea_ln[e]["b"]
                if np.isfinite(m) and np.isfinite(b) and len(inv_T) >= 2:
                    inv_T_fit = np.linspace(inv_T.min(), inv_T.max(), 50)
                    ax.plot(inv_T_fit, np.exp(m * inv_T_fit + b),
                            linestyle="--", linewidth=0.7,
                            color=colors[i], alpha=0.7)

        ax.set_yscale("log")
        ax.set_title(r"$\ln\tau$ fit (inv\_sig)")
        ax.set_xlabel(r"$1/T$")
        ax.set_ylabel(r"$\langle \tau (e_{IS},T) \rangle$")
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")

        norm = Normalize(vmin=min(energies_sorted), vmax=max(energies_sorted))
        sm   = ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(sm, cax=cax).set_label(r"$e_{IS}$")

        # ── Panel 2: total counts per temperature ─────────────────────────────────
        ax = axes[1]
        total_per_temp = []
        for t in all_temps:
            total = 0
            for e in energies_sorted:
                data = tau_data[e]
                if t in data["temps"]:
                    idx = data["temps"].index(t)
                    c   = data["counts"][idx]
                    if np.isfinite(c):
                        total += int(c)
            total_per_temp.append(total)

        ax.bar(range(len(all_temps)), total_per_temp, alpha=0.8)
        if min_total > 0:
            ax.axhline(min_total, color="red", ls="--", lw=0.8,
                       label=f"min\\_total={min_total}")
            ax.legend(fontsize=5)
        ax.set_xticks(range(len(all_temps)))
        ax.set_xticklabels([f"{t:.3f}" for t in all_temps],
                           rotation=45, ha="right", fontsize=5)
        ax.set_yscale("log")
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$\sum N$")
        ax.set_title(r"Total counts per $T$")
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")

        # ── Panel 3: heatmap counts(e_IS, T) ─────────────────────────────────────
        ax = axes[2]
        count_matrix = np.zeros((len(energies_sorted), len(all_temps)))
        for i, e in enumerate(energies_sorted):
            data = tau_data[e]
            for j, t in enumerate(all_temps):
                if t in data["temps"]:
                    idx = data["temps"].index(t)
                    c   = data["counts"][idx]
                    count_matrix[i, j] = int(c) if np.isfinite(c) else 0

        vmin = max(1, count_matrix[count_matrix > 0].min()) \
               if np.any(count_matrix > 0) else 1
        im = ax.imshow(count_matrix, aspect="auto", origin="lower",
                       cmap="viridis",
                       norm=LogNorm(vmin=vmin, vmax=count_matrix.max()),
                       extent=[-0.5, len(all_temps) - 0.5,
                               -0.5, len(energies_sorted) - 0.5])

        ax.set_xticks(range(len(all_temps)))
        ax.set_xticklabels([f"{t:.3f}" for t in all_temps],
                           rotation=45, ha="right", fontsize=5)
        ax.set_yticks(range(len(energies_sorted)))
        ax.set_yticklabels([f"{e:.4f}" for e in energies_sorted], fontsize=5)
        ax.set_xlabel(r"$T$")
        ax.set_ylabel(r"$e_{IS}$")
        ax.set_title(r"Counts per $(e_{IS}, T)$")

        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax).set_label(r"$N$")

        plt.tight_layout()
        LifetimesPlotter._save(fig, save)
        plt.show(block=False)
        plt.pause(20)
        plt.close()

    @staticmethod
    def check_plot_activ(activ_data:dict,save:str="check_activ",energies:list=None):
        """
        
        """
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize

        plt = LifetimesPlotter._setup_plotting()
        available = sorted(activ_data.keys())

        if energies is not None:
            energies_plot = sorted(set(
                min(available,key=lambda x: abs(x-e)) for e in energies
            ))
        else:
            energies_plot = available
        
        weight_names = ["sqrt_n","inv_sig","uniform"]
        weight_labels = {
            "sqrt_n": r"$w = \sqrt{N}$",
            "inv_sig": r"$w = 1/\sigma(\ln\tau)$",
            "uniform": r"$w = 1$",
        }

        fig1,axes1 = plt.subplots(3,3,figsize=(3*3.39,6.0))

        for col,wname in enumerate(weight_names):
            e_arr = []
            m_arr = []
            m_err = []
            R2_arr = []
            chi2_arr = []

            for e in energies_plot:
                entry = activ_data[e]
                if entry is None: continue
                f = entry["fits"].get(wname)
                if f is None: continue
                e_arr.append(e)
                m_arr.append(f["m"])
                m_err.append(f["m_err"])
                R2_arr.append(f["R2"])
                chi2_arr.append(f["chi2_red"])
            
            e_arr = np.array(e_arr)
            m_arr = np.array(m_arr)
            m_err = np.array(m_err)
            R2_arr = np.array(R2_arr)
            chi2_arr = np.array(chi2_arr)

            valid = np.isfinite(m_arr) & np.isfinite(R2_arr)

            ax = axes1[0,col]
            if np.any(valid):
                norm = Normalize(vmin=np.nanmin(R2_arr[valid]),
                                 vmax=np.nanmax(R2_arr[valid]))
                sc = ax.scatter(e_arr[valid],m_arr[valid],c=R2_arr[valid],
                                cmap="viridis",norm=norm,s=15,zorder=3)
                ax.errorbar(e_arr[valid],m_arr[valid],yerr=m_err[valid],
                            fmt="none",ecolor="gray",capsize=2,capthick=0.7,lw=0.7,zorder=2)
                fig1.colorbar(sc,ax=ax).set_label(r"$R^2$")
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$E_a$")
            ax.set_title(weight_labels[wname])
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")

            ax = axes1[1,col]
            if np.any(valid):
                ax.plot(e_arr[valid],R2_arr[valid],marker="o",markersize=3,linestyle="-",
                        label=r"$R^2$")
                ax2 = ax.twinx()
                ax2.plot(e_arr[valid],chi2_arr[valid],marker="s",markersize=3,
                         linestyle="--", color="orange",label=r"$\chi^2_\mathrm{red}$")
                ax2.set_ylabel(r"$\chi^2_\mathrm{red}$")
                ax2.tick_params(labelsize=5)
            ax.axhline(1.0,color="green",lw=0.7,ls="--",alpha=0.6)
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$R^2$")
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")

            ax = axes1[2,col]
            if np.any(valid):
                n_pts_arr = np.array([
                    activ_data[e]["fits"][wname]["n_pts"]
                    if (activ_data[e] is not None and
                        activ_data[e]["fits"].get(wname) is not None)
                        else 0 
                        for e in e_arr[valid] 
                ])
                ax.bar(e_arr[valid],n_pts_arr,
                       width=(e_arr[1]-e_arr[0]) * 0.8 if len(e_arr) > 1 else 0.001,
                       alpha=0.8)
                ax.axhline(2,color="red",ls="--",lw=0.8)
            ax.set_xlabel(r"$e_{IS}$")
            ax.set_ylabel(r"$n_T$")
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")
            ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True,nbins=3))

        plt.tight_layout()
        LifetimesPlotter._save(fig1,save + "_overview")
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        fig2,axes2 = plt.subplots(4,3,figsize=(3*3.39,8.0))

        for col,wname in enumerate(weight_names):
            ranked = []
            for e in energies_plot:
                entry = activ_data[e]
                if entry is None: continue
                f = entry["fits"].get(wname)
                if f is None or not np.isfinite(f["R2"]): continue
                ranked.append((f["R2"],e,entry,f))
            if not ranked: continue
            ranked.sort(key=lambda x: x[0])
            best = ranked[-1]
            worst = ranked[0]

            for row_offset, (R2,e,entry,f) in enumerate([best,worst]):
                tag = "best" if row_offset == 0 else "worst"
                row_f = row_offset * 2
                row_r = row_offset * 2 + 1

                inv_T = entry["inv_T"]
                ln_tau = entry["ln_tau"]
                counts = entry["counts"]

                c_norm = counts / counts.max() if counts.max() > 0 else counts
                scatter_colors = plt.cm.viridis(c_norm)

                ax = axes2[row_f,col]
                ax.scatter(inv_T,ln_tau,c=scatter_colors,s=20,zorder=3)
                ax.plot(f["inv_T_fine"],f["ln_tau_fine"],color="red",lw=1.2,zorder=2)
                ax.set_xlabel(r"$1/T$")
                ax.set_ylabel(r"$\ln\tau$")
                ax.set_title(
                    rf"{weight_labels[wname]} - {tag}"
                    f"\n$e_{{IS}}={e:.5f}$, "
                    rf"$R^2={R2:.3f}$, "
                    rf"$\chi^2_{{\mathrm{{red}}}}={f['chi2_red']:.2f}$"
                )
                ax.xaxis.set_ticks_position("bottom")
                ax.yaxis.set_ticks_position("left")

                ax = axes2[row_r,col]
                ax.axhline(0,color="gray",lw=0.8)
                ax.scatter(inv_T,f["residuals"],c=scatter_colors,s=15,zorder=2)
                if wname == "inv_sig": 
                    sigma_band = 1.0 / f["weights"]
                    ax.fill_between(inv_T,-sigma_band,sigma_band,alpha=0.2,color="red")
                ax.set_xlabel(r"$1/T$")
                ax.set_ylabel(r"$\ln\tau - \ln\tau_\mathrm{fit}$")
                ax.xaxis.set_ticks_position("bottom")
                ax.yaxis.set_ticks_position("left")

            plt.tight_layout()
        LifetimesPlotter._save(fig2,save + "_detail")
        plt.show(block=False)
        plt.pause(5)
        plt.close()
        return

    @staticmethod
    def plot_quality(results:dict,save:str="quality"):
        """
        
        """
        plt = LifetimesPlotter._setup_plotting()

        min_counts = sorted(set(k[0] for k in results.keys()))
        min_totals = sorted(set(k[1] for k in results.keys()))

        metrics = {
            "mean_R2": (r"$\langle R^2 \rangle$", "viridis", False),
            "mean_rel_error": (r"$\langle |\sigma_{E_a}/E_a| \rangle$", "magma_r", False),
            "n_valid": (r"$N_\mathrm{valid}$", "cividis", False),
            "mean_chi2": (r"$\langle \chi^2_\mathrm{red} \rangle$", "inferno_r", True),
        }

        fig, axes = plt.subplots(2,2,figsize=(7,6))

        for ax, (key,(label,cmap,use_log)) in zip(axes.flat, metrics.items()):
            matrix = np.full((len(min_counts),len(min_totals)), np.nan)

            for i, mc in enumerate(min_counts):
                for j, mt in enumerate(min_totals):
                    if (mc,mt) in results:
                        matrix[i,j] = results[(mc,mt)][key]

            if use_log:
                from matplotlib.colors import LogNorm
                valid_vals = matrix[np.isfinite(matrix) & (matrix > 0)]
                if len(valid_vals) > 0:
                    norm = LogNorm(vmin=valid_vals.min(), vmax=valid_vals.max())
                else:
                    norm = None
                im = ax.imshow(matrix,aspect="auto",origin="lower",
                               cmap=cmap,norm=norm)
            else:
                im = ax.imshow(matrix,aspect="auto",origin="lower",
                               cmap=cmap)
                
            ax.set_xticks(range(len(min_totals)))
            ax.set_xticklabels(min_totals)
            ax.set_yticks(range(len(min_counts)))
            ax.set_yticklabels(min_counts)
            ax.set_xlabel("min\\_total")
            ax.set_ylabel("min\\_count")
            ax.set_title(label)

            for i in range(len(min_counts)):
                for j in range(len(min_totals)):
                    val = matrix[i,j]
                    if np.isfinite(val):
                        fmt = f"{val:.3f}" if val < 100 else f"{val:.0f}"
                        ax.text(j,i,fmt,ha="center",va="center",fontsize=5,color="white")
            
            fig.colorbar(im,ax=ax,shrink=0.8)
        plt.tight_layout()
        LifetimesPlotter._save(fig,save)
        plt.show(block=False)
        plt.pause(3)
        plt.close()
        return

# ── IBMAnalysis ───────────────────────────────────────────────────────────────

class IBMAnalysis:
    """
    Interval Bisection Method on a LAMMPS IS-energy trajectory.
    """
    def __init__(self, path):
        self.times, self.e_is = (np.atleast_1d(x) for x in np.loadtxt(path, unpack=True))
        sort         = np.argsort(self.times)
        self.times   = self.times[sort]
        self.e_is    = self.e_is[sort]
        self.e_is    = np.round(self.e_is, 16)
        self.path    = path

        _, unique_idx = np.unique(self.times[::-1], return_index=True)
        unique_idx    = len(self.times) - 1 - unique_idx
        unique_idx.sort()
        self.times = self.times[unique_idx]
        self.e_is  = self.e_is[unique_idx]

    def calc_jumps(self):
        """
        Prints the timepoints between IS-energy jumps.
        """
        for i, eis in enumerate(self.e_is):
            if i == 0:
                continue
            if (self.e_is[i-1] != eis) and (abs(self.times[i-1] - self.times[i]) != 1):
                print(f"{int(self.times[i-1])} {int(self.times[i])}")

    def calc_jumps_short(self):
        """
        
        """
        for i in range(1,len(self.e_is)):
            if self.e_is[i-1] == self.e_is[i]: continue
            if abs(self.times[i-1] - self.times[i]) == 1: continue
            e_before = self.e_is[i-1]
            if e_before in self.e_is[i+1:]: continue
            print(f"{int(self.times[i-1])} {int(self.times[i])}")

    def clean_file(self):
        """
        Removes unnecessary timepoints from the IS trajectory file.
        """
        keep = set()
        keep.add(0)
        keep.add(len(self.e_is) - 1)
        for i, eis in enumerate(self.e_is):
            if i == 0 or i == (len(self.e_is) - 1):
                continue
            if self.e_is[i-1] != eis:
                keep.add(i)
                keep.add(i-1)

        keep        = sorted(keep)
        times_clean = self.times[keep]
        eis_clean   = self.e_is[keep]

        with open(self.path, "w") as f:
            for t, e in zip(times_clean, eis_clean):
                f.write(f"{int(t)} {e:.15f}\n")


# ── MBAnalysis ────────────────────────────────────────────────────────────────

class MBAnalysis:
    """
    Meta-Basin Algorithm on a LAMMPS IS-energy trajectory.
    """
    def __init__(self, path):
        self.times, self.e_is = (np.atleast_1d(x) for x in np.loadtxt(path, unpack=True))
        sort       = np.argsort(self.times)
        self.times = self.times[sort]
        self.e_is  = self.e_is[sort]
        self.e_is  = np.round(self.e_is, 10)
        self.first  = None
        self.last   = None
        self.uniq   = None
        self.nr_int = 0

    @classmethod
    def read_file(cls, path):
        """
        Reads an already computed MB result file (first last e_IS).

        Parameters
        ----------
        path : str
        """
        data = np.loadtxt(path)
        if data.ndim == 1:
            data = data[None, :]
        mb         = cls.__new__(cls)
        mb.times   = None
        mb.e_is    = None
        mb.first   = data[:, 0]
        mb.last    = data[:, 1]
        mb.uniq    = data[:, 2]
        mb.nr_int  = len(mb.uniq)
        return mb
    
    @staticmethod
    def calc_intervals_fast(mb):
        """
        Performs the calculation of intervals faster than before.
        """
        t0 = time.time()
        print("calc_intervals: sortiere und gruppiere...", flush=True)

        e = mb.e_is
        t = mb.times

        sort_e   = np.argsort(e, kind='stable')
        e_sorted = e[sort_e]
        t_sorted = t[sort_e]

        # Grenzen zwischen Gruppen gleicher Energie
        boundaries = np.flatnonzero(np.diff(e_sorted)) + 1
        starts = np.concatenate([[0], boundaries])
        ends   = np.concatenate([boundaries, [len(e_sorted)]])

        # first/last = min/max der Zeit innerhalb jeder Gruppe
        # Nutzt reduceat statt Python-Loop
        mb.first  = np.minimum.reduceat(t_sorted, starts)
        mb.last   = np.maximum.reduceat(t_sorted, starts)
        mb.uniq   = e_sorted[starts]
        mb.nr_int = len(starts)

        mb._sort_intervals()
        print(f"calc_intervals: {mb.nr_int} Intervalle in {time.time()-t0:.1f}s", flush=True)

    @staticmethod
    def mb_rehwald_timed(mb):
        stages = [
            ("calc_intervals",    lambda: MBAnalysis.calc_intervals_fast(mb)),
            ("calc_cut (1/2)",    lambda: mb.calc_cut_rehwald()),
            ("calc_cut (2/2)",    lambda: mb.calc_cut_rehwald()),
            ("calc_combine",      lambda: mb.calc_combine()),
            ("calc_remove",       lambda: mb.calc_remove()),
            ("check",             lambda: (mb._check_duplicates(), mb._check_overlaps())),
        ]
        total_t = time.time()
        bar = tqdm(stages, desc="mb_rehwald", unit="phase")
        for name, fn in bar:
            bar.set_postfix(phase=name, nr_int=getattr(mb, 'nr_int', '?'))
            t0 = time.time()
            fn()
            dt = time.time() - t0
            tqdm.write(f"  ✓ {name:<22} {mb.nr_int:>8} Intervalle   {dt:6.1f}s")
        tqdm.write(f"  Gesamt: {time.time()-total_t:.1f}s")

    def mb_rehwald(self):
        """
        MB algorithm according to C. Rehwald (2012).
        """
        self.calc_intervals()
        self.calc_cut_rehwald()
        self.calc_cut_rehwald()
        self.calc_combine()
        self.calc_remove()
        self._check_duplicates()
        self._check_overlaps()

    def calc_intervals(self):
        """
        Calculates intervals from first to last occurrence of each minimum.
        """
        uniq       = np.unique(self.e_is)
        first_list = []
        last_list  = []
        uniq_list  = []
        for eis in uniq:
            idx = np.where(self.e_is == eis)[0]
            first_list.append(self.times[idx[0]])
            last_list.append(self.times[idx[-1]])
            uniq_list.append(eis)
        self.first  = np.array(first_list)
        self.last   = np.array(last_list)
        self.uniq   = np.array(uniq_list)
        self.nr_int = len(uniq_list)
        self._sort_intervals()

    def calc_cut_rehwald(self):
        """
        Cuts intervals overlapping less than 50% (Rehwald 2012).
        """
        updates = {}
        for i in range(self.nr_int):
            for j in range(i+1, self.nr_int):
                if self.first[j] >= self.last[i]: break
                if self.first[i] < self.first[j] < self.last[i] < self.last[j]:
                    overlap = self._calc_overlap(i, j)
                    if overlap <= 0.5:
                        updates[j] = max(updates.get(j, self.first[j]), self.last[i])
        for j, new_first in updates.items():
            self.first[j] = new_first
        self._sort_intervals()

    def calc_combine(self):
        """
        Combines intervals with an overlap greater than 50%.
        """
        parent = list(range(self.nr_int))

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py
                return True
            return False

        for i in range(self.nr_int):
            for j in range(i+1, self.nr_int):
                if self.first[j] >= self.last[i]: break
                if self.first[i] < self.first[j] < self.last[i] < self.last[j]:
                    if self._calc_overlap(i, j) > 0.5:
                        union(i, j)

        clusters = {}
        for i in range(self.nr_int):
            clusters.setdefault(find(i), []).append(i)

        new_first, new_last, new_uniq = [], [], []
        for members in clusters.values():
            if len(members) == 1:
                idx = members[0]
                new_first.append(self.first[idx])
                new_last.append(self.last[idx])
                new_uniq.append(self.uniq[idx])
            else:
                new_first.append(np.min(self.first[members]))
                new_last.append(np.max(self.last[members]))
                new_uniq.append(np.min(self.uniq[members]))

        self.first  = np.array(new_first)
        self.last   = np.array(new_last)
        self.uniq   = np.array(new_uniq)
        self.nr_int = len(new_first)
        self._sort_intervals()

    def calc_remove(self):
        """
        Removes intervals that lie completely within other intervals.
        """
        remove = set()
        for i in range(self.nr_int):
            if i in remove:
                continue
            for j in range(self.nr_int):
                if i == j or j in remove:
                    continue
                if self.first[i] <= self.first[j] and self.last[i] >= self.last[j]:
                    self.uniq[i] = min(self.uniq[i], self.uniq[j])
                    remove.add(j)

        keep_mask       = np.ones(self.nr_int, dtype=bool)
        keep_mask[list(remove)] = False
        self.first  = self.first[keep_mask]
        self.last   = self.last[keep_mask]
        self.uniq   = self.uniq[keep_mask]
        self.nr_int = int(np.sum(keep_mask))
        self._sort_intervals()

    def _calc_overlap(self, i, j):
        max_len = max(self.last[i] - self.first[i], self.last[j] - self.first[j])
        overlap = max(self.first[i], self.first[j]) - min(self.last[i], self.last[j])
        return overlap / max_len

    def _sort_intervals(self):
        sort       = np.argsort(self.first)
        self.first = self.first[sort]
        self.last  = self.last[sort]
        self.uniq  = self.uniq[sort]
        return

    def _check_duplicates(self):
        unique_energies = set(self.uniq)
        if len(unique_energies) < self.nr_int:
            counts     = {}
            for e in self.uniq:
                counts[e] = counts.get(e, 0) + 1
            duplicates = {e: c for e, c in counts.items() if c > 1}
            print(f"WARNING: {len(duplicates)} energies appear multiple times.")
        return
    
    def _check_overlaps(self):
        overlaps = []
        for i in range(self.nr_int):
            for j in range(i+1, self.nr_int):
                left  = max(self.first[i], self.first[j])
                right = min(self.last[i],  self.last[j])
                if right > left:
                    overlaps.append((i, j, right - left))
        if overlaps:
            print(f"WARNING: {len(overlaps)} intervals overlap.")
        return
    

    def write_results(self,path:str):
        """
        Writes the resulting intervals to a file.

        Parameters
        ----------
        path : str
        """
        with open(path, "w") as f:
            for i, eis in enumerate(self.uniq):
                f.write(f"{self.first[i]} {self.last[i]} {eis}\n")
        return
    
    @staticmethod
    def plot_MBs(path:str,save:str):
        """
        
        """
        from matplotlib.collections import LineCollection
        import matplotlib as mpl
        mpl.use("TkAgg")
        import scienceplots
        plt.style.use(["science","nature"])
        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "sans-serif"

        mb = MBAnalysis.read_file(path=path)
        fig,ax = plt.subplots(figsize=(6,3))
        
        ax.hlines(mb.uniq,mb.first,mb.last)
        ax.set_xlabel(r"MD-Steps")
        ax.set_ylabel(r"$e_{\mathrm{MB}}$")
        plt.tight_layout()
        fig.savefig(f"{save}.png",dpi=300)
        plt.show(block=False)
        plt.pause(3)
        plt.close()
        return

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    """
    :Author: Simon Georg Kellers
    :Email: s_kell14@uni-muenster.de
    :Date: 2026
    :Version: 6.0.0
    """
    parser = argparse.ArgumentParser(
        description=(
            "Performs different analyses on minimized IS-energies:\n"
            "  Lifetimes   : Calculates MB lifetimes from MB result files.\n"
            "  Fitting     : Fits tau(e_IS) from a single Lifetimes result file.\n"
            "  TauT        : Gathers tau(e_IS,T) from multiple Lifetimes result files.\n"
            "  Life_Plot   : Plots tau(e_IS) and phi(e_IS) from a Lifetimes result file.\n"
            "  TauT_Plot   : Plots tau(e_IS,T) from a TauT result file.\n"
            "  Activ_Plot  : Plots E_a(e_IS) from a TauT result file.\n"
            "  Global_Plot : Plots global mean tau vs T from Lifetimes result files.\n"
            "  IBM         : Detects timepoints between IS-energy changes.\n"
            "  MB          : Identifies meta-basins from IS trajectory.\n"
            "  MB_Plot     : Plots the MB lifetimes against the trajectory.\n"
            "  Clean       : Removes unnecessary timepoints from IS trajectory.\n"
            "  Debug_Plot  : Plots tau(e_IS) and phi(e_IS) from a Lifetimes result file.\n"
            "  Dist_Plot   : Plots the number of data per T and e_IS.\n"
            "  Counts      : Writes Counts per Energy Bin from a TauT file.\n"
            "  Quality     : Evaluates the Quality of the E_A Fit based on min-count and min-total.\n"
        ),
        epilog="For more details see the documentation.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("-i", "--input",   type=str, required=True,
                        help="Path to input file (or glob pattern).")
    parser.add_argument("-o", "--output",  type=str, default=None,
                        help="Path to output file.")
    parser.add_argument("-m", "--mode",    type=str, required=True,
                        choices=["Lifetimes", "Fitting", "TauT", "Counts", "Quality", "Activ",
                                 "Life_Plot", "TauT_Plot", "Activ_Plot", "Global_Plot",
                                 "Debug_Plot", "Dist_Plot", "MB_Plot",
                                 "IBM", "MB", "Clean",
                                 "Check_Binned", "Check_Fitted", "Check_TauT", "Check_Activ",
                                 "MB_Timed", "IBM_Short"])
    parser.add_argument("-e", "--energies", nargs="+", type=float, default=None,
                        help="Energies at which to evaluate fits.")
    parser.add_argument("-t", "--thresh",   type=float, default=1.0,
                        help="Threshold multiplier for tau filtering.")
    parser.add_argument("-b","--min-count",type=int,default=1,
                        help="Minimum occupation for a bin to be considered.")
    parser.add_argument("-B","--min-total",type=int,default=0,
                        help="Min total counts across all T.")

    args = parser.parse_args()

    # ── Lifetimes ─────────────────────────────────────────────────────────────
    if args.mode == "Lifetimes":
        files = [f for f in glob.glob(args.input) if f.endswith(".dat")]
        if not files:
            print(f"ERROR: No files found for pattern '{args.input}'")
            sys.exit(1)

        tau_dict = {}

        for f in files:
            mb = MBAnalysis.read_file(f)
            for i, eis in enumerate(mb.uniq):
                tau = float(mb.last[i] - mb.first[i])
                if eis not in tau_dict:
                    tau_dict[eis] = []
                tau_dict[eis].append(tau)

        eis_sorted = sorted(tau_dict.keys())
        life = Lifetimes(
            eis      = eis_sorted,
            tau_mean = [float(np.mean(tau_dict[e])) for e in eis_sorted],
            tau_std  = [float(np.std(tau_dict[e]))  for e in eis_sorted],
            dist     = [len(tau_dict[e])             for e in eis_sorted],
        )
        life.calc_mean()
        LifetimesIO.write_results(life, path=args.output)

    # ── Fitting ───────────────────────────────────────────────────────────────
    elif args.mode == "Fitting":
        life = LifetimesIO.read_file(args.input)
        life.filter_tau(thresh=args.thresh)

        fitter = LifetimesFitter(life)
        fitter.fit_ln_tau()
        fitter.fit_dir_tau(energies=args.energies,
                           min_count=args.min_count)

        if args.output:
            LifetimesIO.write_fits(fitter, path=args.output, energies=args.energies,
                                   min_count=args.min_count)

    # ── TauT ──────────────────────────────────────────────────────────────────
    elif args.mode == "TauT":
        tau_data = LifetimesIO.gather_tau_vs_T(
            path=args.input, energies=args.energies, thresh=args.thresh,
            min_count=args.min_count
        )
        if args.output:
            LifetimesIO.write_tau_vs_T(tau_data, path=args.output)

    # ── Life_Plot ─────────────────────────────────────────────────────────────
    elif args.mode == "Life_Plot":
        life = LifetimesIO.read_file(args.input)
        life.filter_tau(thresh=args.thresh)

        fitter = LifetimesFitter(life)
        fitter.fit_ln_tau()
        fitter.fit_dir_tau(energies=args.energies,
                           min_count=args.min_count)
        fitter.binned_tau(energies=args.energies,
                          min_count=args.min_count)

        plotter = LifetimesPlotter(life, fitter)
        plotter.plot_dist(save=args.output or "dist")

    # ── TauT_Plot ─────────────────────────────────────────────────────────────
    elif args.mode == "TauT_Plot":
        tau_data = LifetimesIO.read_tau_vs_T(args.input)
        ea_ln  = LifetimesFitter.calc_activation(tau_data, mode="Fit",
                                                 min_total=args.min_total)
        ea_dir = LifetimesFitter.calc_activation(tau_data, mode="Dir",
                                                 min_total=args.min_total)
        ea_bin = LifetimesFitter.calc_activation(tau_data, mode="Bin",
                                                 min_total=args.min_total)
        LifetimesPlotter.plot_tau_vs_T(tau_data, save=args.output or "tau_vs_T",
                                       energies=args.energies,min_total=args.min_total,
                                       ea_ln=ea_ln,ea_dir=ea_dir,ea_bin=ea_bin)

    # ── Activ_Plot ────────────────────────────────────────────────────────────
    elif args.mode == "Activ_Plot":
        tau_data = LifetimesIO.read_tau_vs_T(args.input)

        ea_ln  = LifetimesFitter.calc_activation(tau_data, mode="Fit",
                                                 min_total=args.min_total)
        ea_dir = LifetimesFitter.calc_activation(tau_data, mode="Dir",
                                                 min_total=args.min_total)
        ea_bin = LifetimesFitter.calc_activation(tau_data, mode="Bin",
                                                 min_total=args.min_total)

        LifetimesPlotter.plot_activation(
            ea_ln=ea_ln, ea_dir=ea_dir, ea_bin=ea_bin, save=args.output or "activation",
            energies=args.energies, min_total=args.min_total,tau_data=tau_data
        )

    # ── Global_Plot ───────────────────────────────────────────────────────────
    elif args.mode == "Global_Plot":
        LifetimesPlotter.plot_global_tau(
            path=args.input, save=args.output or "global_tau"
        )

    # ── IBM ───────────────────────────────────────────────────────────────────
    elif args.mode == "IBM":
        mb = IBMAnalysis(path=args.input)
        mb.calc_jumps()

    elif args.mode == "IBM_Short":
        mb = IBMAnalysis(path=args.input)
        mb.calc_jumps_short()

    # ── MB ────────────────────────────────────────────────────────────────────
    elif args.mode == "MB":
        mb = MBAnalysis(path=args.input)
        mb.mb_rehwald()
        mb.write_results(path=args.output)

    elif args.mode == "MB_Timed":
        mb = MBAnalysis(path=args.input)
        mb.mb_rehwald_timed(mb=mb)
        mb.write_results(path=args.output)

    # ── Clean ─────────────────────────────────────────────────────────────────
    elif args.mode == "Clean":
        mb = IBMAnalysis(path=args.input)
        mb.clean_file()

    # ── MB_Plot ─────────────────────────────────────────────────────────────────
    elif args.mode == "MB_Plot":
        MBAnalysis.plot_MBs(path=args.input,save=args.output)

    # ── Debug_Plot ────────────────────────────────────────────────────────────
    elif args.mode == "Debug_Plot":
        life = LifetimesIO.read_file(args.input)
        life.filter_tau(thresh=args.thresh)

        fitter = LifetimesFitter(life)
        fitter.fit_ln_tau()
        fitter.fit_dir_tau(energies=args.energies,
                           min_count=args.min_count)
        fitter.binned_tau(energies=args.energies,
                          min_count=args.min_count)

        LifetimesPlotter.plot_debug(life, fitter, save=args.output or "debug")

    # ── Dist_Plot ─────────────────────────────────────────────────────────────
    elif args.mode == "Dist_Plot":
        tau_data = LifetimesIO.read_tau_vs_T(args.input)
        LifetimesPlotter.plot_counts(tau_data=tau_data,save=args.output or "counts",
                                     energies=args.energies,min_total=args.min_total)
        
    # ── Counts ────────────────────────────────────────────────────────────────
    elif args.mode == "Counts":
        tau_data = LifetimesIO.read_tau_vs_T(args.input)
        LifetimesIO.write_counts(
            tau_data=tau_data,path=args.output or "counts.dat",
            energies=args.energies,min_total=args.min_total
        )

    # ── Quality ────────────────────────────────────────────────────────────────
    elif args.mode == "Quality":
        results = LifetimesFitter.calc_quality(
            path=args.input,
            energies=args.energies,
            thresh=args.thresh,
            mode="Dir",
        )
        LifetimesPlotter.plot_quality(results,save=args.output or "quality")

    elif args.mode == "Check_Binned":
        life = LifetimesIO.read_file(args.input)
        LifetimesPlotter.check_plot_binned(
            life,
            min_count=args.min_count, save=args.output
        )

    elif args.mode == "Check_Fitted":
        life = LifetimesIO.read_file(args.input)
        LifetimesPlotter.check_plot_fitted(
            life,
            min_count=args.min_count,save=args.output
        )

    elif args.mode == "Check_TauT":
        tau_data = LifetimesIO.check_tauT(
            path=args.input,
            energies=args.energies,
            min_count=args.min_count,
        )
        LifetimesPlotter.check_plot_tauT(
            tau_data=tau_data,
            save=args.output or "check_tauT",
            energies=args.energies,
            min_total=args.min_total,
        )

    elif args.mode == "Check_Activ":
        tau_data = LifetimesIO.check_tauT(
            path=args.input,
            energies=args.energies,
            min_count=args.min_count,
        )
        activ_data = LifetimesFitter.check_activ(
            tau_data=tau_data,
            min_total=args.min_total,
        )
        LifetimesPlotter.check_plot_activ(
            activ_data=activ_data,
            save=args.output or "check_activ",
            energies=args.energies,
        )

    elif args.mode == "Activ":
        tau_data = LifetimesIO.check_tauT(
            path=args.input, 
            energies=args.energies,
            min_count=args.min_count,
        )
        activ_data = LifetimesFitter.check_activ(
            tau_data=tau_data,
            min_total=args.min_total,
        )
        LifetimesIO.write_activ(
            activ_data=activ_data,
            path=args.output or "activ.dat",
            energies=args.energies,
        )

if __name__ == "__main__":
    main()
