"""Arrow and Polars helpers for Bamboo."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pyarrow as pa


def to_polars(table: pa.Table) -> Any:
    """Convert a PyArrow table to a Polars DataFrame."""
    import polars as pl

    return pl.from_arrow(table)