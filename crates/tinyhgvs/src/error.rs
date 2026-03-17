//! Structured parser errors returned by [`crate::parse_hgvs`].
//!
//! The parser separates failure into a small number of user-facing categories:
//!
//! - invalid syntax that does not match the supported subset at all
//! - unsupported syntax that is biologically recognizable HGVS, but not yet
//!   modeled by `tinyhgvs`
//! - semantic constraints reserved for future higher-level validation
//!
//! In practice, the most useful fields are usually:
//!
//! - [`ParseHgvsError::code`] for a stable machine-readable diagnostic such as
//!   `unsupported.allele`
//! - [`ParseHgvsError::message`] for a short explanation
//! - [`ParseHgvsError::fragment`] for the most relevant unsupported fragment
//! - [`ParseHgvsError::parser_version`] for tracing the crate release that
//!   produced the error

use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// A structured error returned when an HGVS string cannot be parsed.
///
/// The error describes both what failed and how `tinyhgvs` classified the
/// failure.
///
/// This is especially useful for syntaxes that are valid HGVS but not yet
/// supported by the current data model. For example, an allele expression such
/// as `NC_000001.11:g.[123G>A;345del]` returns
/// `code == "unsupported.allele"` rather than a generic parse failure.
///
/// # Examples
///
/// Invalid syntax:
///
/// ```rust
/// use tinyhgvs::{ParseHgvsErrorKind, parse_hgvs};
///
/// let error = parse_hgvs("not-an-hgvs-string").unwrap_err();
/// assert_eq!(error.kind(), ParseHgvsErrorKind::InvalidSyntax);
/// assert_eq!(error.code(), "invalid.syntax");
/// ```
///
/// Unsupported but recognized HGVS syntax:
///
/// ```rust
/// use tinyhgvs::{ParseHgvsErrorKind, parse_hgvs};
///
/// let error = parse_hgvs("p.(Gln576SerfsTer21)").unwrap_err();
/// assert_eq!(error.kind(), ParseHgvsErrorKind::UnsupportedSyntax);
/// assert_eq!(error.code(), "unsupported.protein_frameshift");
/// assert_eq!(error.fragment(), Some("fs"));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseHgvsError {
    kind: ParseHgvsErrorKind,
    code: &'static str,
    message: String,
    input: String,
    fragment: Option<String>,
    parser_version: &'static str,
}

/// High-level classes of parse failures exposed by `tinyhgvs`.
///
/// # Examples
///
/// ```rust
/// use tinyhgvs::{ParseHgvsErrorKind, parse_hgvs};
///
/// let invalid = parse_hgvs("bad").unwrap_err();
/// assert_eq!(invalid.kind(), ParseHgvsErrorKind::InvalidSyntax);
///
/// let unsupported = parse_hgvs("NM_004006.3:r.spl").unwrap_err();
/// assert_eq!(unsupported.kind(), ParseHgvsErrorKind::UnsupportedSyntax);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParseHgvsErrorKind {
    /// The input does not match the grammar supported by the parser.
    InvalidSyntax,
    /// The input matches a known HGVS family that is intentionally unsupported.
    UnsupportedSyntax,
    /// The input parsed structurally but violates a semantic constraint.
    SemanticConstraint,
}

impl ParseHgvsError {
    pub(crate) fn new(
        kind: ParseHgvsErrorKind,
        code: &'static str,
        message: impl Into<String>,
        input: impl Into<String>,
        fragment: Option<String>,
    ) -> Self {
        Self {
            kind,
            code,
            message: message.into(),
            input: input.into(),
            fragment,
            parser_version: env!("CARGO_PKG_VERSION"),
        }
    }

    pub(crate) fn invalid(input: &str) -> Self {
        Self::new(
            ParseHgvsErrorKind::InvalidSyntax,
            "invalid.syntax",
            "failed to parse HGVS variant",
            input,
            None,
        )
    }

    pub(crate) fn unsupported(
        code: &'static str,
        message: &'static str,
        input: &str,
        fragment: Option<String>,
    ) -> Self {
        Self::new(
            ParseHgvsErrorKind::UnsupportedSyntax,
            code,
            message,
            input,
            fragment,
        )
    }

    /// Returns the broad error class.
    pub fn kind(&self) -> ParseHgvsErrorKind {
        self.kind
    }

    /// Returns the machine-friendly diagnostic code such as `unsupported.allele`.
    pub fn code(&self) -> &'static str {
        self.code
    }

    /// Returns the human-readable error message.
    pub fn message(&self) -> &str {
        &self.message
    }

    /// Returns the full HGVS string that failed to parse.
    pub fn input(&self) -> &str {
        &self.input
    }

    /// Returns the most relevant fragment recognized by the diagnostic layer.
    ///
    /// For example, a protein frameshift may return `"fs"`, while an allele
    /// expression may return a bracketed fragment such as `"[123G>A;345del]"`.
    pub fn fragment(&self) -> Option<&str> {
        self.fragment.as_deref()
    }

    /// Returns the `tinyhgvs` crate version that produced the error.
    pub fn parser_version(&self) -> &'static str {
        self.parser_version
    }
}

impl Display for ParseHgvsError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}] {}: `{}` (tinyhgvs {})",
            self.code, self.message, self.input, self.parser_version
        )
    }
}

impl Error for ParseHgvsError {}
