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
from operator import attrgetter
from warnings import warn
import re
compartment_finder = re.compile(r"^\s*(\[[A-Za-z]\])\s*:*")
# Regular expressions to match the arrows
_reversible_arrow_finder = re.compile("<(-+|=+)>")
_forward_arrow_finder = re.compile("(-+|=+)>")
_reverse_arrow_finder = re.compile("<(-+|=+)")

class Reaction():
    def __init__(
        self,
        id: str = None,
        name: str = "",
        subsystem: str = "",
        metabolites: Dict[str, float] = {},
        gene_reaction_rule: str = "",
        lower_bound: float = 0.0,
        upper_bound: float = 1000.0,
        kcat: float = 100000,
        mw: float = 1.0,
        kcat_MW: float = None,
        g0: float = -99999.0,
        annotation: dict = {}
    ):
        """Initialize a simple object with an identifier.

        Parameters
        ----------
        id: string, optional
            the identifier to associate with the object
        name: string, optional
            The name to associate with the object. Default "".
        subsystem: string, optional
            The subsystem to associate with the object. Default "".
        lower_bound: float, optional
            The lower bound to associate with the object. Default 0.0.
        upper_bound: float, optional
            The upper bound to associate with the object. Default 0.0.
        kcat: float, optional
            The kcat to associate with the object. Default 100000.
        mw: float, optional
            The molecular weight to associate with the object. Default 1.0.
        g0: float, optional
            The g0 to associate with the object. Default -99999.0.
        **kwargs: Any
            Keyword arguments passed on to `json.dumps`.

        Objects will have notes and _annotation as dicitionaries, initialized as empty
        dictionaries.
        """
        self.id = id
        self.name = name
        self.subsystem = subsystem
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.kcat = kcat
        self.mw = mw
        self.kcat_MW = kcat_MW
        self.g0 = g0
        self.gene_reaction_rule = gene_reaction_rule
        self.metabolites = metabolites
        self.notes = {}
        self.annotation = annotation
    
    def __str__(self):
        return f"{self.id}: {self.build_reaction_string(model,use_metabolite_names=False)}"

    def build_reaction_string(self,model,use_metabolite_names: bool = False) -> str:
        """Generate a human readable reaction str.

        Parameters
        ----------
        use_metabolite_names: bool
            Whether to use metabolite names (when True) or metabolite ids (when False,
            default).

        Returns
        -------
        str
            A human readable str.
        """

        def format(number: Union[int, float]) -> str:
            return "" if number == 1 else str(number).rstrip(".") + " "

        reactant_bits = []
        product_bits = []
        for key,value in self.metabolites.items():
            if value >= 0:
                product_bits.append(format(value) + key if not use_metabolite_names else model.metabolites[key].name)
            else:
                reactant_bits.append(format(abs(value)) + key if not use_metabolite_names else model.metabolites[key].name)

        reaction_string = " + ".join(reactant_bits)
        if not self.reversibility:
            if self.lower_bound < 0 and self.upper_bound <= 0:
                reaction_string += " <-- "
            else:
                reaction_string += " --> "
        else:
            reaction_string += " <=> "
        reaction_string += " + ".join(product_bits)
        return reaction_string
    
    def set_id(self, value: str) -> None:
        """Set the the standard gibbs free energy.

        Parameters
        ----------
        value: float
            The value to set the the standard gibbs free energy.

        Setting the the standard gibbs free energy (float) for the construction of thermodynamics constraint model.
        """
        # Validate bounds before setting them.
        self.id = value

    @property
    def get_kcat(self):
        """Get the kcat.

        Returns
        -------
        float
            The kcat of the reaction.
        """
        return self.kcat

    def set_kcat(self, value: float) -> None:
        """Set the kcat.

        Parameters
        ----------
        value: float
            The value to set the kcat.

        Setting the kcat for the construction of enzyme constraint model.
        """
        # Validate bounds before setting them.
        if value <= 0:
            raise ValueError(
            f"kcat can not < 0"
            )
        self.kcat = value

    @property
    def get_mw(self) -> float:
        """Get the molecule weight.

        Returns
        -------
        float
            The molecule weight of the reaction.
        """
        return self.mw

    def set_mw(self, value: float) -> None:
        """Set the molecule weight.

        Parameters
        ----------
        value: float
            The value to set the molecule weight.

        Setting the molecule weight (float) for the construction of enzyme constraint model.
        """
        # Validate bounds before setting them.
        if value <= 0:
            raise ValueError(
            f"MW can not < 0"
            )
        self.mw = value

    @property
    def get_g0(self) -> float:
        """Get the standard gibbs free energy.

        Returns
        -------
        float
            The the standard gibbs free energy of the reaction.
        """
        return self.g0
    
    def set_g0(self, value: float) -> None:
        """Set the the standard gibbs free energy.

        Parameters
        ----------
        value: float
            The value to set the the standard gibbs free energy.

        Setting the the standard gibbs free energy (float) for the construction of thermodynamics constraint model.
        """
        # Validate bounds before setting them.
        self.g0 = value

    def get_upper_bound(self) -> float:
        """Get the upper bound.

        Returns
        -------
        float
            The upper bound of the reaction.
        """
        return self.upper_bound
    
    def get_lower_bound(self) -> float:
        """Get the lower bound.

        Returns
        -------
        float
            The lower bound of the reaction.
        """
        return self.lower_bound
    
    def set_upper_bound(self, value: float) -> None:
        """Set the upper bound.

        Parameters
        ----------
        value: float
            The value to set the upper bound.

        Setting the upper bound (float) for the construction of enzyme constraint model.
        """
        # Validate bounds before setting them.
        self.check_bounds(self.lower_bound, value)
        self.upper_bound = value
    
    def set_lower_bound(self, value: float) -> None:
        """Set the lower bound.

        Parameters
        ----------
        value: float
            The value to set the lower bound.

        Setting the lower bound (float) for the construction of enzyme constraint model.
        """
        # Validate bounds before setting them.
        self.check_bounds(value, self.upper_bound)
        self.lower_bound = value

    @property
    def bounds(self) -> Tuple[float, float]:
        """Get or the bounds.

        Returns
        -------
        tuple: lower_bound, upper_bound
            A tuple of floats, representing the lower and upper bound.
        """
        return self.lower_bound, self.upper_bound

    @bounds.setter
    def bounds(self, value: Union[Tuple[float, float], Sequence[float]]) -> None:
        """Set the bounds directly, using a tuple or list.

        Parameters
        ----------
        value: tuple or sequence
            The lower bound and upper bound. Invalid bounds will raise ValueError.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Raises
        ------
        ValueError
            If lower bound higher than upper bound, via _check_bounds.

        """
        lower, upper = value
        # Validate bounds before setting them.
        self.check_bounds(lower, upper)
        self.lower_bound = lower
        self.upper_bound = upper

    @property
    def reactants(self) -> List[str]:
        """Return a list of reactants for the reaction.

        Returns
        -------
        list
            A list of the metabolites consudmed (coefficient < 0) by the reaction.
        """
        return [k for k, v in self.metabolites.items() if v < 0]

    @property
    def products(self) -> List[str]:
        """Return a list of products for the reaction.

        Returns
        -------
        list
            A list of the metabolites produced (coefficient > 0) by the reaction.
        """
        return [k for k, v in self.metabolites.items() if v >= 0]

    def check_bounds(self,lb: float, ub: float) -> None:
        """Check if the lower and upper bounds are valid.

        Arguments
        ----------
        lb: float
            The lower bound.
        ub: float
            The upper bound.

        Raises
        ------
        ValueError
            If the lower bound is higher than upper bound.
        """
        if lb > ub:
            raise ValueError(
                f"The lower bound must be less than or equal to the upper bound "
                f"({lb} <= {ub})."
            )    
        
    def reversibility(self) -> bool:
        """Whether the reaction can proceed in both directions (reversible).

        This is computed from the current upper and lower bounds.

        Returns
        -------
        bool:
            True if the reaction is reversible (lower bound lower than 0 and upper
            bound is higher than 0).

        """
        return self.lower_bound < 0 < self.upper_bound
    
    def add_metabolites(self, reaction_dict):
        pop_list = []
        for key, value in reaction_dict.items():
            if key in self.metabolites.keys():
                self.metabolites[key] += value
            else:
                self.metabolites[key] = value
        for key, value in self.metabolites.items():
            if value == 0:
                pop_list.append(key)
        for key in pop_list:
            self.metabolites.pop(key)

    @property
    def genes(self):
        gpr = self.gene_reaction_rule
        genes_list =[]
        for gene in gpr.split(" "):
            if gene != "and":
                genes_list.append(gene)
        return genes_list

    def get_coefficient(self, metabolite_id):
        if metabolite_id in self.metabolites.keys():
            return self.metabolites[metabolite_id]
        else:
            raise ValueError(f"{metabolite_id} not in {self.id}")

    def build_reaction_from_string(
        self,
        reaction_str: str,
        verbose: bool = True,
        fwd_arrow: Optional[AnyStr] = None,
        rev_arrow: Optional[AnyStr] = None,
        reversible_arrow: Optional[AnyStr] = None,
        term_split: str = "+",
    ) -> None:
        """Build reaction from reaction equation reaction_str using parser.

        Takes a string and using the specifications supplied in the optional
        arguments infers a set of metabolites, metabolite compartments and
        stoichiometries for the reaction.  It also infers the reversibility
        of the reaction from the reaction arrow.

        Changes to the associated model are reverted upon exit when using
        the model as a context.

        Parameters
        ----------
        reaction_str : str
            a string containing a reaction formula (equation)
        verbose: bool
            setting verbosity of function (default True)
        fwd_arrow : AnyStr, optional
            Str or bytes that encode forward irreversible reaction arrows (default
            None).
        rev_arrow : AnyStr, optional
            Str or bytes that encode backward irreversible reaction arrows (default
            None).
        reversible_arrow : AnyStr, optional
            Str or bytes that encode reversible reaction arrows (default None).
        term_split : str
            dividing individual metabolite entries (default "+")".

        Raises
        ------
        ValueError
            No arrow found in reaction string.
        """
        # set the arrows
        forward_arrow_finder = (
            _forward_arrow_finder
            if fwd_arrow is None
            else re.compile(re.escape(fwd_arrow))
        )
        reverse_arrow_finder = (
            _reverse_arrow_finder
            if rev_arrow is None
            else re.compile(re.escape(rev_arrow))
        )
        reversible_arrow_finder = (
            _reversible_arrow_finder
            if reversible_arrow is None
            else re.compile(re.escape(reversible_arrow))
        )
        found_compartments = compartment_finder.findall(reaction_str)
        if len(found_compartments) == 1:
            compartment = found_compartments[0]
            reaction_str = compartment_finder.sub("", reaction_str)
        else:
            compartment = ""

        # reversible case
        arrow_match = reversible_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            self.bounds = (0,1000)
        else:  # irreversible
            # try forward
            arrow_match = forward_arrow_finder.search(reaction_str)
            if arrow_match is not None:
                self.bounds = (0,1000)
            else:
                # must be reverse
                arrow_match = reverse_arrow_finder.search(reaction_str)
                if arrow_match is None:
                    raise ValueError(f"no suitable arrow found in '{reaction_str}'")
                else:
                    self.bounds = (-1000, 0)
        reactant_str = reaction_str[: arrow_match.start()].strip()
        product_str = reaction_str[arrow_match.end() :].strip()

        # self.subtract_metabolites(self.metabolites, combine=True)

        for substr, factor in ((reactant_str, -1), (product_str, 1)):
            if len(substr) == 0:
                continue
            for term in substr.split(term_split):
                term = term.strip()
                if term.lower() == "nothing":
                    continue
                if " " in term:
                    num_str, met_id = term.split()
                    num = float(num_str.lstrip("(").rstrip(")")) * factor
                else:
                    met_id = term
                    num = factor
                met_id += compartment
                # try:
                #     met = model.metabolites.get_by_id(met_id)
                # except KeyError:
                #     if verbose:
                #         print(f"unknown metabolite '{met_id}' created")
                #     met = Metabolite(met_id)
                self.add_metabolites({met_id: num})
    # @property
    # def objective_coefficient(self) -> float:
    #     if self.id == c
