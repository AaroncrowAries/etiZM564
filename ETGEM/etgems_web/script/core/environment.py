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
class Environment():
    def __init__(
        self,
        p_h: float = 7.5,
        p_mg: float = 3.0,
        ionic_strength: str = "0.1M",
        temperature: str = "310.15K",
    ):
        self.p_h = p_h
        self.p_mg = p_mg
        self.ionic_strength = ionic_strength
        self.temperature = temperature
    
    def get_p_h(self):
        """Get the pH."""
        return self.p_h
    
    def get_p_mg(self):
        """Get the pMg."""
        return self.p_mg
    
    def get_ionic_strength(self):
        """Get the ionic strength."""
        return self.ionic_strength
    
    def get_temperature(self):
        """Get the temperature."""
        return self.temperature
    
    def set_p_h(self, p_h: float) -> None:
        """Set the pH.

        Arguments
        ----------
        p_h: float
            The pH.
        """
        self.p_h = p_h
    
    def set_p_mg(self, p_mg: float) -> None:
        """Set the pMg.

        Arguments
        ----------
        p_mg: float
            The pMg.
        """
        self.p_mg = p_mg
    
    def set_ionic_strength(self, ionic_strength: str) -> None:
        """Set the ionic strength.

        Arguments
        ----------
        ionic_strength: str
            The ionic strength.
        """
        self.ionic_strength = ionic_strength

    def set_temperature(self, temperature: str) -> None:
        """Set the temperature.

        Arguments
        ----------
        temperature: str
            The temperature.
        """
        self.temperature = temperature
    
    