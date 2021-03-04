class Error(Exception):
    """Base class for exceptions in this project."""

class DtaFileNotFound(Error):
    """Raised when a DTASelect file is not found at the expected path."""

class CimageFlatFileNotFound(Error):
    """Raised when a cimage flatfile (output_to_excel.txt) is not found at the expected path."""

class FilterNotFoundException(Error):
    """Raised when filter is not found."""
