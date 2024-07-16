"""
Functions and classes used by other scripts

Author: Diana Szeliova
Lat update: 15.7.2024
"""

import cobra
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

random.seed(42)


class Model:
    """
    class for manipulation of cobra models
    """

    def __init__(self, cobra_model):
        self.model = cobra_model
        self.model.objective = "0"
        self.reset_bounds()
        self.rate_ids = {
            "67": "q_MSG",
            "96": "q_Glc",
            "132": "q_Trehalose",
            "47": "q_Lactate",
            "106": "q_EtOH",
            "104": "q_Citrate",
            "54": "q_Alanine",
            "60": "q_Arginine",
            "62": "q_Asparagine",
            "63": "q_Aspartate",
            "69": "q_Histidine",
            "71": "q_Isoleucine",
            "73": "q_Leucine",
            "75": "q_Methionine",
            "77": "q_Phenylalanine",
            "79": "q_Proline",
            "80": "q_Serine",
            "82": "q_Threonine",
            "84": "q_Valine",
            "108": "q_Glycine",
            "992": "q_Lysine",
        }

    def set_lower_bound(self, rxn, bound):
        """set reaction lower bound"""
        self.model.reactions.get_by_id(rxn).lower_bound = bound

    def set_upper_bound(self, rxn, bound):
        """set reaction upper bound"""
        self.model.reactions.get_by_id(rxn).upper_bound = bound

    def set_bounds(self, rxn, lower_bound, upper_bound):
        """set lower and upper bound of reaction"""
        self.model.reactions.get_by_id(rxn).lower_bound = lower_bound
        self.model.reactions.get_by_id(rxn).upper_bound = upper_bound

    def reset_bounds(self):
        """
        Reset all bounds in the model do defaults (0, 1000) for irreversible,
        (-1000, 1000) for reversible reactions.
        """
        for rxn in self.model.reactions:
            rxn.upper_bound = 1000
            rxn.lower_bound = [0, -1000][rxn.reversibility]

    def identify_exchanges(self):
        """
        Identify exchange reactions -- reactions that have no reactants or
        products.
        Returns:
            list of exchange reaction IDs
        """
        exchanges = []
        for r in self.model.reactions:
            if not r.products or not r.reactants:
                exchanges.append(r.id)

        return exchanges

    def identify_transports(self, reaction_ids):
        """
        Identify transport reactions based on their names.
        Arguments:
            reaction_ids: a dictionary of reaction ID: reaction name
        Returns:
            list of transporter IDs
        """
        transports = []
        for r in self.model.reactions:
            reaction_name = reaction_ids[f"R_{r.id}"]
            if ("trans_" in reaction_name) or ("transporter" in reaction_name):
                transports.append(r.id)

        return transports

    def make_transporters_reversible(self, reaction_ids):
        """
        Make all transporters in the model reversible.
        """
        transports = self.identify_transports(reaction_ids)
        for transport_id in transports:
            self.set_lower_bound(transport_id, -1000)
            self.set_upper_bound(transport_id, 1000)

    def print_reaction_with_names(self, reaction_id, metabolite_ids, reaction_ids):
        """
        Print reaction string with metabolite names.
        Arguments:
            reaction_id: ID of reaction to be printed
            metabolite_ids: dictionary of metabolite ID: metabolite name
            reaction_ids: a dictionary of reaction ID: reaction name
        """
        rxn = self.model.reactions.get_by_id(reaction_id)
        reaction_string = str(rxn)
        split_rxn = reaction_string.split(" ")

        for part in split_rxn:
            try:
                name = metabolite_ids[f"M_{part}"]
                reaction_string = reaction_string.replace(f" {part}", f" {name}")
            except KeyError:
                pass

        reaction_name = reaction_ids[f"R_{reaction_id}"]
        reaction_string = reaction_string.replace(f"{reaction_id}:", f"{reaction_id} ({reaction_name}):\n")

        print(reaction_string)

    def set_rates(self, rates_mean, rates_sd, mode="mean"):
        """
        Set the lower and upper bounds of reactions according to the provided
        data and the chosen mode.
        Arguments:
            rates_mean: data frame with reaction rates
            rates_mean: data frame with rate SDs
            mode: mean - takes mean exchange rate as lower and upper bounds
                mean_sd - mean minus SD as lower bound, mean plus SD as upper bound
                sampled - random sample from [rate-2*SD, rate+2*SD]
        """

        for rxn_id, rxn_name in self.rate_ids.items():
            try:
                mean_rate = rates_mean[rxn_name]
                sd_rate = rates_sd[rxn_name]
            except KeyError:
                continue

            if mode == "mean":
                lb = mean_rate
                ub = mean_rate
            if mode == "mean_sd":
                lb = mean_rate-sd_rate
                ub = mean_rate+sd_rate
            if mode == "sampled":
                bound1 = random.uniform(mean_rate-2*sd_rate, mean_rate+2*sd_rate)
                bound2 = random.uniform(mean_rate-2*sd_rate, mean_rate+2*sd_rate)
                lb, ub = sorted([bound1, bound2])

            self.set_lower_bound(rxn_id, lb)
            self.set_upper_bound(rxn_id, ub)

    def run_pfba_sampled(self, n_feasible, fluxes, growth_rates):
        """
        Run pFBA and add results to existing objects
        Arguments:
            n_feasible: integer with the number of feasible solutions so far
            fluxes: list of solution objects from pFBA
            growth rates: list of growth rates
        Returns:
            a tuple of updated n_feasible, fluxes and growth rates
        """
        mu = self.model.slim_optimize()
        if mu == mu:
            n_feasible += 1
            fba = cobra.flux_analysis.pfba(self.model)
            fluxes.append(fba.fluxes)
            growth_rates.append(fba.fluxes[0])

        return (n_feasible, fluxes, growth_rates)


def get_central_fluxes(central_metabolism, fba_solution):
    """returns a list of fluxes for reactions in 'central_metabolism'
    Arguments:
        central_metabolism: data_frame with reaction IDs that we are interested in.
            If one reaction is mapped to several IDs, they are
            separated by a comma
        fba_solution: data frame with predicted fluxes in the column "fluxes"
    Returns:
        mapped_fluxes: list of mapped fluxes
    """

    mapped_fluxes = []

    for index, row in central_metabolism.iterrows():
        if row.ID != row.ID:
            mapped_fluxes.append(0)
            continue

        rxn_ids = str(row.ID).split(", ")
        parallel_fluxes = []
        for rxn_id in rxn_ids:
            parallel_fluxes.append(fba_solution.fluxes[rxn_id])

        # if multiple IDs map to the same reaction, sum up the fluxes
        mapped_fluxes.append(sum(parallel_fluxes))

    return mapped_fluxes


def plot_sampling(fba_dictionary, experimental_mean, experimental_sd, feature,
                  ylabel=None):
    """
    Plot results from sampled pFBA for one reaction as boxplot
    + experimental data as a point with error bars
    Arguments:
        fba_dictionary: dictionary {condition: list of fluxes of target reaction}
        experimental_mean: pd.DataFrame with mean experimental values
        experimental_sd: pd.DataFrame with SD of experimental values
        feature: feature to plot
        ylabel: optional, label for the y-axis. If not provided, feature is used
    """
    fig, ax = plt.subplots()
    ax.boxplot(fba_dictionary.values())
    ticks = [name.replace(" Cell Density", "") for name in fba_dictionary.keys()]
    ax.set_xticklabels(ticks)
    sim = mpatches.Patch(color='darkorchid', label='Simulation')

    # sampled data
    for i, key in enumerate(fba_dictionary):
        y = fba_dictionary[key]
        x = np.random.normal(1+i, 0.06, size=len(y))
        plt.plot(x, y, color="darkorchid", marker='.', 
                 alpha=0.2, linestyle="")

    ys = [i+1.3 for i in range(len(fba_dictionary))]

    # experimental data
    mean = experimental_mean[feature][fba_dictionary.keys()]
    sd2 = experimental_sd[feature][fba_dictionary.keys()]*2
    plot2 = ax.scatter(ys, mean,
                       color="mediumseagreen",
                       label="Experiment")
    plt.errorbar(ys, mean, sd2,
                 linestyle="",
                 color="mediumseagreen")

    plt.ylim(0, max(mean)+max(sd2)*1.5)
    plt.ylabel(feature)
    if ylabel:
        plt.ylabel(ylabel)
    plt.legend(handles=[sim, plot2], loc="upper left")
