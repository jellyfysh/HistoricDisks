# HistoricDisks - Synopsis of pressure data, sampling algorithms and pressure estimators for the hard-disk model of
# statistical physics
# https://github.com/jellyfysh/HistoricDisks
# Copyright (C) 2022 Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, and Werner Krauth
#
# This file is part of HistoricDisks.
#
# HistoricDisks is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# HistoricDisks is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with HistoricDisks in the LICENSE file.
# If not, see <https://www.gnu.org/licenses/>.
#
# If you use HistoricDisks in published work, please cite the following reference (see [Li2022] in References.bib):
# Botao Li, Yoshihiko Nishikawa, Philipp Höllmer, Louis Carillo, A. C. Maggs, Werner Krauth;
# Hard-disk pressure computations—a historic perspective.
# J. Chem. Phys. 21 December 2022; 157 (23): 234111. https://doi.org/10.1063/5.0126437
#
"""Executable Python script that parses csv files and converts between different units."""
import csv
from enum import Enum
import math
from typing import List


class XUnit(Enum):
    """
    Enumeration class for the possible units on the x-axis.
    """

    SPECIFIC_VOLUME = "V / V_0"  # String "V / V_0" can be accessed via XUnit.SPECIFIC_VOLUME.value.
    VOLUME_PER_PARTICLE = "v / (2 sigma)^2"
    PACKING_FRACTION = "packing fraction"
    REDUCED_DENSITY = "reduced density"

    @staticmethod
    def from_string(value: str) -> 'XUnit':
        """
        Return the enumeration member with the given value.

        Parameters
        ----------
        value : str
            The value.

        Returns
        -------
        XUnit
            The enumeration member.

        Raises
        ------
        RuntimeError
            If no enumeration member has the given value.
        """
        for x_unit in XUnit:
            if x_unit.value == value:
                return x_unit
        raise RuntimeError(f"Unknown unit on the x-axis: {value}")


class YUnit(Enum):
    """
    Enumeration class for the possible units on the y-axis.
    """
    LEGACY = "beta P V_0 / N"
    MODERN = "beta P (2 sigma)^2"

    @staticmethod
    def from_string(value: str) -> 'YUnit':
        """
        Return the enumeration member with the given value.

        Parameters
        ----------
        value : str
            The value.

        Returns
        -------
        YUnit
            The enumeration member.

        Raises
        ------
        RuntimeError
            If no enumeration member has the given value.
        """
        for y_unit in YUnit:
            if y_unit.value == value:
                return y_unit
        raise RuntimeError(f"Unknown unit on the y-axis: {value}")


class Tag(object):
    """
    Class to store the tag of data.

    Attributes
    ----------
    number_disks : str
        The number of disks.
    aspect_ratio : str
        The aspect ratio.
    x_unit : XUnit
        The unit on the x-axis.
    y_unit : YUnit
        The unit on the y-axis.
    type : str
        The type.
    """
    def __init__(self, number_disks: str, aspect_ratio: str, x_unit: XUnit, y_unit: YUnit, type: str) -> None:
        """
        Initialize the tag.

        Parameters
        ----------
        number_disks : str
            The number of disks.
        aspect_ratio : str
            The aspect ratio.
        x_unit : XUnit
            The unit on the x-axis.
        y_unit : YUnit
            The unit on the y-axis.
        type : str
            The type.
        """
        self.number_disks = number_disks
        self.aspect_ratio = aspect_ratio
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.type = type

    def __str__(self) -> str:
        """
        Return the string representation of the tag (e.g., for printing).

        Returns
        -------
        str
            The string representation.
        """
        return (f"Tag(number_disks={self.number_disks}, aspect_ratio={self.aspect_ratio}, x_unit={self.x_unit}, " +
                f"y_unit={self.y_unit}, type={self.type})")


class Data(object):
    """
    Class to store data.

    Attributes
    ----------
    tag : Tag
        The tag of the data.
    x_values : typing.List[float]
        The values on the x-axis.
    y_values : typing.List[float]
        The values on the y-axis.
    y_errors : typing.List[float]
        The errors of the values on the y-axis.
    """
    def __init__(self, tag: Tag, x_values: List[float], y_values: List[float],
                 y_errors: List[float]) -> None:
        """
        Initialize the data.

        Attributes
        ----------
        tag : Tag
            The tag of the data.
        x_values : typing.List[float]
            The values on the x-axis.
        y_values : typing.List[float]
            The lower values on the y-axis.
        y_errors : typing.List[float]
            The upper values on the y-axis.
        """
        self.tag = tag
        self.x_values = x_values
        self.y_values = y_values
        self.y_errors = y_errors

    def __str__(self) -> str:
        """
        Return the string representation of the data (e.g., for printing).

        Returns
        -------
        str
            The string representation.
        """
        print_string = f"Tag: {self.tag}"
        print_string += f"\nx values: {self.x_values}"
        print_string += f"\ny values: {self.y_values}"
        print_string += f"\ny errors: {self.y_errors}"
        return print_string


class DataGroup(object):
    """
    Class to store a data group from a parsed csv file.

    A data group consists in a title of the csv file, and several data instances.

    Attributes
    ----------
    title : str
        The title.
    data : List[Data]
        The data.
    """
    def __init__(self, title: str, data: List[Data]) -> None:
        """
        Initialize the data group.

        Parameters
        ----------
        title : str
            The title.
        data : List[Data]
            The data.
        """
        self.title = title
        self.data = data

    def __str__(self) -> str:
        """
        Return the string representation of the data group (e.g., for printing).

        Returns
        -------
        str
            The string representation.
        """
        print_string = "------------" + "-" * len(self.title)
        print_string += f"\nData group {self.title}."
        for data in self.data:
            print_string += f"\n{data}"
        print_string += "\n------------" + "-" * len(self.title)
        return print_string


def read(filename: str, title: str = "") -> DataGroup:
    """
    Parse a data group from the given csv file.

    If no title is given as a parameter to this function, it tries to parse it from the csv file.

    Each data group is stored in a csv file. The header of the csv file contains information regarding the source of the
    data. This function reads the value of the field "label" as the default title of the data group and ignore all other
    fields in the header.

    Each data group is composed by several data subgroups, specified by the setup of computation (number of disks,
    aspect ratio, etc.) and the representation of the data (unit of volumes, unit of pressures, etc.). The subgroups are
    separated by an empty line. Each subgroup has a three-line header providing information regarding the simulation
    from which the data is collected. This function reads in the fields:
        N: number of disks
        Aspect ratio:the aspect ratio of the simulation box
        X unit: the unit for x-coordinate (either volume or density)
        Y unit: the unit for pressure
        type: indicating whether the error of the pressure is given by the size of the error bar or originated from the
        marker size in the figure
    Other fields are not read into the program. For a detailed description of the format of the csv files, please check
    the Readme of the digitized data.

    Parameters
    ----------
    filename : str
        The filename of the csv file.
    title : str, optional
        The title of the data group.

    Returns
    -------
    DataGroup
        The populated data group with the data from the csv file.

    Raises
    ------
    AssertionError
        If the specified format of the csv file is not met.
    """
    data = []
    label = None
    with open(filename, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',', strict=True)
        for row in csv_reader:
            if row[0].startswith("Label = "):
                assert not row[1] and not row[2]  # Assert that empty.
                assert label is None  # Assert that label only set once.
                label = row[0].partition("Label = ")[2]
                assert label  # Assert that not empty.
                continue

            # Start of new data.
            if row[0].startswith("N = "):
                number_disks = row[0].partition("N = ")[2]
                assert number_disks
                aspect_ratio = row[1].partition("Aspect ratio = ")[2]
                assert aspect_ratio
                row = next(csv_reader)  # Get next row.
                x_unit_value = row[0].partition("X unit = ")[2]
                assert x_unit_value
                x_unit = XUnit.from_string(x_unit_value)
                y_unit_value = row[1].partition("Y unit = ")[2]
                assert y_unit_value
                y_unit = YUnit.from_string(y_unit_value)
                type = row[2].partition("type = ")[2]
                assert type
                # Parse data until empty line.
                x_values = []
                lower_y_values = []
                upper_y_values = []
                # We can use the same csv_reader iterator to get the following rows in the csv file.
                for sub_row in csv_reader:
                    # Empty row.
                    if not sub_row[0] and not sub_row[1] and not sub_row[2]:
                        break
                    x_values.append(float(sub_row[0]))
                    lower_y_values.append(float(sub_row[1]))
                    upper_y_values.append(float(sub_row[2]))
                tag = Tag(number_disks, aspect_ratio, x_unit, y_unit, type)
                data.append(Data(tag, x_values, lower_y_values, upper_y_values))
    if title:
        label = title
    assert label
    data_group = DataGroup(label, data)
    return data_group


def convert_x(x_values: List[float], actual_unit: XUnit, wanted_unit: XUnit) -> List[float]:
    """
    Convert the given values on the x-axis in the given unit, and return a list of values in the new unit.

    Parameters
    ----------
    x_values : List[float]
        The values on the x-axis.
    actual_unit : XUnit
        The unit of the given values.
    wanted_unit : XUnit
        The unit of the returned values.

    Returns
    -------
    List[float]
        The values on the x-axis in the new unit.

    Raises
    ------
    AssertionError
        If one of the units is unknown.
    """
    if actual_unit is XUnit.SPECIFIC_VOLUME:
        if wanted_unit is XUnit.SPECIFIC_VOLUME:
            conversion_factor = 1.0
            invert = False
        elif wanted_unit is XUnit.VOLUME_PER_PARTICLE:
            conversion_factor = math.sqrt(3.0) / 2.0
            invert = False
        elif wanted_unit is XUnit.PACKING_FRACTION:
            conversion_factor = math.pi * math.sqrt(3.0) / 6.0
            invert = True
        else:
            assert wanted_unit is XUnit.REDUCED_DENSITY
            conversion_factor = 2.0 / math.sqrt(3.0)
            invert = True
    elif actual_unit is XUnit.VOLUME_PER_PARTICLE:
        if wanted_unit is XUnit.SPECIFIC_VOLUME:
            conversion_factor = 2.0 / math.sqrt(3.0)
            invert = False
        elif wanted_unit is XUnit.VOLUME_PER_PARTICLE:
            conversion_factor = 1.0
            invert = False
        elif wanted_unit is XUnit.PACKING_FRACTION:
            conversion_factor = math.pi / 4.0
            invert = True
        else:
            assert wanted_unit is XUnit.REDUCED_DENSITY
            conversion_factor = 1.0
            invert = True
    elif actual_unit is XUnit.PACKING_FRACTION:
        if wanted_unit is XUnit.SPECIFIC_VOLUME:
            conversion_factor = math.pi * math.sqrt(3.0) / 6.0
            invert = True
        elif wanted_unit is XUnit.VOLUME_PER_PARTICLE:
            conversion_factor = math.pi / 4.0
            invert = True
        elif wanted_unit is XUnit.PACKING_FRACTION:
            conversion_factor = 1.0
            invert = False
        else:
            assert wanted_unit is XUnit.REDUCED_DENSITY
            conversion_factor = 4.0 / math.pi
            invert = False
    else:
        assert actual_unit is XUnit.REDUCED_DENSITY
        if wanted_unit is XUnit.SPECIFIC_VOLUME:
            conversion_factor = 2.0 / math.sqrt(3.0)
            invert = True
        elif wanted_unit is XUnit.VOLUME_PER_PARTICLE:
            conversion_factor = 1.0
            invert = True
        elif wanted_unit is XUnit.PACKING_FRACTION:
            conversion_factor = math.pi / 4.0
            invert = False
        else:
            assert wanted_unit is XUnit.REDUCED_DENSITY
            conversion_factor = 1.0
            invert = False
    if invert:
        return [conversion_factor / value for value in x_values]
    else:
        return [conversion_factor * value for value in x_values]


def convert_y(y_values: List[float], actual_unit: YUnit, wanted_unit: YUnit) -> List[float]:
    """
    Convert the given values on the y-axis in the given unit, and return a list of values in the new unit.

    Parameters
    ----------
    y_values : List[float]
        The values on the y-axis.
    actual_unit : XUnit
        The unit of the given values.
    wanted_unit : XUnit
        The unit of the returned values.

    Returns
    -------
    List[float]
        The values on the y-axis in the new unit.

    Raises
    ------
    AssertionError
        If one of the units is unknown.
    """
    if actual_unit is YUnit.LEGACY:
        if wanted_unit is YUnit.LEGACY:
            conversion_factor = 1.0
        else:
            assert wanted_unit is YUnit.MODERN
            conversion_factor = 2.0 / math.sqrt(3.0)
    else:
        assert actual_unit is YUnit.MODERN
        if wanted_unit is YUnit.LEGACY:
            conversion_factor = math.sqrt(3.0) / 2.0
        else:
            assert wanted_unit is YUnit.MODERN
            conversion_factor = 1.0
    return [conversion_factor * value for value in y_values]


def convert_data(data_group: DataGroup, x_unit: XUnit, y_unit: YUnit) -> None:
    """
    Convert the data in the data group to the given units on the x- and y-axis.

    Parameters
    ----------
    data_group : DataGroup
        The data group.
    x_unit : XUnit
        The desired unit on the x-axis.
    y_unit : YUnit
        The desired unit on the y-axis.
    """
    for data in data_group.data:
        data.x_values = convert_x(data.x_values, data.tag.x_unit, x_unit)
        data.tag.x_unit = x_unit
        data.y_values = convert_y(data.y_values, data.tag.y_unit, y_unit)
        data.y_errors = convert_y(data.y_errors, data.tag.y_unit, y_unit)
        data.tag.y_unit = y_unit


def main() -> None:
    """
    Parse data from a csv file into a data group, and print it before and after converting units on both axes.
    """
    data = read("../../DigitizedData/Mak_2006/Mak_2006.csv")
    print(data)
    convert_data(data, XUnit.SPECIFIC_VOLUME, YUnit.MODERN)
    print(data)


if __name__ == "__main__":
    main()
