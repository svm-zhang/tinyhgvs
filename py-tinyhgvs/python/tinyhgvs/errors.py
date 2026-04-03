"""Public exception types raised by :mod:`tinyhgvs`."""

from __future__ import annotations

from enum import Enum


class ParseHgvsErrorKind(str, Enum):
    """Broad categories used by :class:`TinyHGVSError`.

    Attributes:
        INVALID_SYNTAX: The input does not match the currently supported HGVS grammar.
        UNSUPPORTED_SYNTAX: The input belongs to a recognized HGVS family that is not supported yet.
        SEMANTIC_CONSTRAINT: Reserved for future semantic validation errors after structural parsing.

    Examples:
        Invalid syntax that does not match the supported grammar:
        >>> from tinyhgvs import TinyHGVSError, parse_hgvs
        >>> try:
        ...     parse_hgvs("not an hgvs variant")
        ... except TinyHGVSError as error:
        ...     error.kind
        <ParseHgvsErrorKind.INVALID_SYNTAX: 'invalid_syntax'>

        Unsupported syntax that is recognized but not implemented:
        >>> try:
        ...     parse_hgvs("NC_000001.11:g.[123G>A;345del]")
        ... except TinyHGVSError as error:
        ...     error.kind
        <ParseHgvsErrorKind.UNSUPPORTED_SYNTAX: 'unsupported_syntax'>
    """

    #: The input does not match the currently supported parser grammar.
    INVALID_SYNTAX = "invalid_syntax"
    #: The input belongs to a known HGVS family that is not supported yet.
    UNSUPPORTED_SYNTAX = "unsupported_syntax"
    #: Reserved for future semantic validation beyond syntax parsing.
    SEMANTIC_CONSTRAINT = "semantic_constraint"


class TinyHGVSError(ValueError):
    """Structured exception for `tinyhgvs` parse failures.

    Attributes:
        kind: Broad error class derived from the Rust error model.
        code: Stable diagnostic code such as ``unsupported.allele``.
        message: Human-readable explanation of the failure.
        input: Original HGVS string that failed to parse.
        fragment: Relevant unsupported fragment when one was detected.
        parser_version: ``tinyhgvs`` version that produced the error.

    Examples:
        Unsupported allele syntax:
        >>> from tinyhgvs import TinyHGVSError, parse_hgvs
        >>> try:
        ...     parse_hgvs("NC_000001.11:g.[123G>A;345del]")
        ... except TinyHGVSError as error:
        ...     (error.kind.value, error.code, error.fragment)
        ('unsupported_syntax', 'unsupported.allele', '[')

        Invalid syntax:
        >>> try:
        ...     parse_hgvs("NM_004006.2:c.5697delA")
        ... except TinyHGVSError as error:
        ...     (error.kind.value, error.code)
        ('invalid_syntax', 'invalid.syntax')

        Unsupported quantified protein insertion syntax:
        >>> try:
        ...     parse_hgvs("p.Arg78_Gly79insXaa[23]")
        ... except TinyHGVSError as error:
        ...     (error.code, error.fragment)
        ('unsupported.protein_insertion_payload', 'Xaa[...]')
    """

    def __init__(
        self,
        kind: ParseHgvsErrorKind | str,
        code: str,
        message: str,
        input: str,
        fragment: str | None,
        parser_version: str,
    ) -> None:
        """Build a Python exception instance from Rust diagnostic details.

        Args:
            kind: Broad error kind.
            code: Stable diagnostic code.
            message: Human-readable error message.
            input: Full HGVS string that failed.
            fragment: Relevant fragment recognized by the diagnostic layer.
            parser_version: ``tinyhgvs`` version that produced the error.
        """
        self.kind = ParseHgvsErrorKind(kind)
        self.code = code
        self.message = message
        self.input = input
        self.fragment = fragment
        self.parser_version = parser_version

        super().__init__(
            f"[{self.code}] {self.message}: `{self.input}` (tinyhgvs {self.parser_version})"
        )


__all__ = ["ParseHgvsErrorKind", "TinyHGVSError"]
