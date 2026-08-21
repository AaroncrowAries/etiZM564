from typing import (
    TYPE_CHECKING,
    AnyStr,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)
from .reaction import Reaction
from .metabolite import Metabolite
from .environment import Environment
from collections import defaultdict


class Model():
    def __init__(
        self,
        genes: dict = {},
        reactions: dict = {},
        metabolites: dict = {},
        environment: Environment = Environment(),
        obj_name: str = "",
        medium: dict = {},
        K_value: float = 0.0,
        E_total: float = 1.0,
        solver: str = "cplex_direct",
    ):
        self.genes = genes
        self.reactions = reactions
        self.metabolites = metabolites
        self.environment = environment
        self.obj_name = obj_name
        self.medium = medium
        self.K_value = K_value
        self.E_total = E_total
        self.solver = solver
    
    def __str__(self):
        return f"{self.id}: {self.obj_name}"

    def get_solver(self):
        """Get the solver."""
        return self.solver
    
    def set_solver(self, solver: str) -> None:
        """Set the solver.

        Arguments
        ----------
        solver: str
            The solver.
        """
        self.solver = solver

    def get_medium(self):
        """Get the medium."""
        return self.medium
    
    def add_medium(self, reaction_id: str, reaction_value: float) -> None:
        """Add a medium.

        Arguments
        ----------
        medium: dict
            The medium.
        """
        self.medium[reaction_id] = reaction_value
    
    def remove_medium(self, reaction_id: str) -> None:
        """Remove a medium.

        Arguments
        ----------
        medium: dict
            The medium.
        """
        del self.medium[reaction_id]

    def get_K_value(self):
        """Get the K value."""
        return self.K_value
    
    def set_K_value(self, K_value: float) -> None:
        """Set the K value.

        Arguments
        ----------
        K_value: float
            The K value.
        """
        self.K_value = K_value
    
    def get_E_total(self):
        """Get the total enzyme."""
        return self.E_total
    
    def set_E_total(self, E_total: float) -> None:
        """Set the total enzyme.

        Arguments
        ----------
        E_total: float
            The total enzyme.
        """
        self.E_total = E_total
    
    def get_obj_name(self):
        """Get the objective name."""
        return self.obj_name
    
    def set_obj_name(self, obj_name: str) -> None:
        """Set the objective name.

        Arguments
        ----------
        obj_name: str
            The objective name.
        """
        self.obj_name = obj_name

    def add_metabolites(self, metabolite: Metabolite):
        """Add a metabolite.

        Arguments
        ----------
        metabolite: Metabolite
            The metabolite.
        """
        if metabolite.id not in self.metabolites.keys():
            self.metabolites[metabolite.id] = metabolite
        else:
            raise ValueError("The metabolite already exists.")

    def remove_metabolites(self, metabolite_id: str):
        """Remove a metabolite.

        Arguments
        ----------
        metabolite_id: str
            The metabolite id.
        """
        del self.metabolites[metabolite_id]

    def add_reactions(self, reaction: Reaction):
        """Add a reaction.

        Arguments
        ----------
        reaction: Reaction
            The reaction.
        """
        if reaction.id not in self.reactions.keys():
            self.reactions[reaction.id] = reaction
            for metabolite in reaction.metabolites.keys():
                if metabolite not in self.metabolites.keys():
                    new_metabolite = Metabolite()
                    new_metabolite.id = metabolite
                    self.add_metabolites(new_metabolite)
        else:
            raise ValueError("The reaction already exists.")
        
    def remove_reactions(self, reaction_id: str):
        """Remove a reaction.

        Arguments
        ----------
        reaction_id: str
            The reaction id.
        """
        del self.reactions[reaction_id]
    
    def get_reaction_by_id(self, reaction_id):
        """Get a reaction by id."""
        return self.reactions[reaction_id]
    
    def get_metabolite_by_id(self, metabolite_id):
        """Get a metabolite by id."""
        return self.metabolites[metabolite_id]
    
    def get_gene_by_id(self, gene_id):
        """Get a gene by id."""
        return self.genes[gene_id]
    
    def get_reactions_by_geneid(self, gene_id):
        """Get a gene by id."""
        result = []
        # for reaction in self.reactions.values():
        for reaction in self.reactions:
            gene_list = reaction.gene_reaction_rule.split(" ")
            if gene_id in gene_list:
                result.append(reaction.id)
        return result
    
    def change_reaction_id(self, old_id, new_id):
        reaction = self.get_reaction_by_id(old_id)
        reaction.id = new_id
        self.remove_reactions(old_id)
        self.add_reactions(reaction)

    def get_objective_coefficient(self, reaction_id):
        self.get_reaction_by_id(reaction_id)
        if reaction_id == self.obj_name:
            return 1
        else:
            return 0

    def check_mass_balance(self,reaction_id):
        reaction_element_dict = defaultdict(int)
        reaction = self.get_reaction_by_id(reaction_id)
        for metabolite_id, coefficient in reaction.metabolites.items():
            metabolite = self.get_metabolite_by_id(metabolite_id)
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += coefficient * metabolite.charge
            if metabolite.elements is None:
                raise ValueError(f"No elements found in metabolite {metabolite.id}")
            for element, amount in metabolite.elements.items():
                reaction_element_dict[element] += coefficient * amount
        return {k: v for k, v in reaction_element_dict.items() if v != 0}

    def merge_reaction(self, reaction_1, reaction_2):
        """Merge two reactions into one reaction.
        Arguments:
        ----------
        * model: a dictionary of model.
        * reaction_1: The first reaction to be merged.
        * reaction_2: The second reaction to be merged.
        Returns:
        ----------
        * new_reaction: The merged reaction.
        """
        new_reaction = Reaction()
        # for reaction in self.reactions.values():
        for reaction in self.reactions:
            if reaction.id == reaction_1:
                first_reaction = reaction
            elif reaction.id == reaction_2:
                second_reaction = reaction
        merged_dict = first_reaction.metabolites.copy()
        for key, value in second_reaction.metabolites.items():
            if key in merged_dict:
                merged_dict[key] += value
                if abs(merged_dict[key]) == 0:
                    del merged_dict[key]
            else:
                merged_dict[key] = value
        if first_reaction.gene_reaction_rule != second_reaction.gene_reaction_rule:
            new_gene_reaction_rule = first_reaction.gene_reaction_rule + " and " + second_reaction.gene_reaction_rule
        else:
            new_gene_reaction_rule = first_reaction.gene_reaction_rule
        if first_reaction.annotation != second_reaction.annotation:
            new_annotation = first_reaction.annotation.copy()
            try:
                second_reaction.annotation["ec-code"]
            except:
                pass
            else:
                try:
                    first_reaction.annotation["ec-code"]
                except:
                    new_annotation["ec-code"] = second_reaction.annotation["ec-code"]
                else:
                    new_annotation["ec-code"] = [first_reaction.annotation["ec-code"],second_reaction.annotation["ec-code"]]
        else:
            new_annotation = first_reaction.annotation.copy()
        new_reaction.id = first_reaction.id + "+" + second_reaction.id
        new_reaction.name = first_reaction.name + "+" + second_reaction.name
        new_reaction.metabolites = merged_dict
        new_reaction.lower_bound = (first_reaction.lower_bound + second_reaction.lower_bound)/2
        new_reaction.upper_bound = (first_reaction.upper_bound + second_reaction.upper_bound)/2
        new_reaction.gene_reaction_rule = new_gene_reaction_rule
        new_reaction.annotation = new_annotation
        new_reaction.set_g0(first_reaction.g0 + second_reaction.g0)
        new_reaction.set_mw(first_reaction.mw + second_reaction.mw)
        new_reaction.set_kcat = min(first_reaction.kcat,second_reaction.kcat)
        self.get_reaction_by_id(reaction_1).bounds = (0,0)
        self.get_reaction_by_id(reaction_2).bounds = (0,0)
        self.add_reactions(new_reaction)
    
    def split_reaction(self, merge_reaction_id):
        merge_reaction = self.get_reaction_by_id(merge_reaction_id)
        merge_reaction.bounds = (0,0)
        [reaction1_id,reaction2_id] = merge_reaction_id.split("+")
        reaction1 = self.get_reaction_by_id(reaction1_id)
        reaction2 = self.get_reaction_by_id(reaction2_id)
        reaction1.bounds = (0,1000)
        reaction2.bounds = (0,1000)

    def add_demand_reaction(self, metabolite_id):
        new_demand_Reaction = Reaction()
        new_demand_Reaction.set_id("DM_"+metabolite_id)
        new_metabolites = dict({metabolite_id: -1})
        new_demand_Reaction.metabolites = new_metabolites
        self.add_reactions(new_demand_Reaction)
    
    def add_sink_reaction(self, metabolite_id):
        new_demand_Reaction = Reaction()
        new_demand_Reaction.set_id("SK_"+metabolite_id)
        new_metabolites = dict({metabolite_id: 1})
        new_demand_Reaction.metabolites = new_metabolites
        self.add_reactions(new_demand_Reaction)