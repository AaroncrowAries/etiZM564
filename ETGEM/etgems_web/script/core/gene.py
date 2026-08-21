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
class Gene():
    def __init__(
        self,
        id: str = None,
        name: str = "",
        uniprot: str = "",
        annotation: dict = {}
    ):
        self.id = id
        self.name = name
        self.uniprot = uniprot
        self.annotation = annotation
    def __repr__(self):
        return f"Gene(id={self.id}, name={self.name}, uniprot={self.uniprot})"
    def __str__(self):
        return f"Gene(id={self.id}, name={self.name}, uniprot={self.uniprot})"
