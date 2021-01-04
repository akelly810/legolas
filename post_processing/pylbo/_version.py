from packaging import version

# The current pylbo version number.
__version__ = "0.4.1"


class LegolasVersion:
    """
    Parser for the Legolas version number, allows for direct comparison

    Parameters
    ----------
    version_number : str
        The Legolas version number in the form "x.x.x".
    """
    def __init__(self, version_number):
        self._validate_version(version_number)
        self._version_number = version.parse(version_number)

    @staticmethod
    def _validate_version(version_number):
        if not isinstance(version_number, str):
            raise ValueError("Version number should be a string.")

    def __lt__(self, other):
        self._validate_version(other)
        return self._version_number < version.parse(other)

    def __le__(self, other):
        self._validate_version(other)
        return self._version_number <= version.parse(other)

    def __gt__(self, other):
        self._validate_version(other)
        return self._version_number > version.parse(other)

    def __ge__(self, other):
        self._validate_version(other)
        return self._version_number >= version.parse(other)

    def __eq__(self, other):
        self._validate_version(other)
        return self._version_number == version.parse(other)

    def __ne__(self, other):
        self._validate_version(other)
        return self._version_number != version.parse(other)
