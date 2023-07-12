## Introduction 
This directory contains the digitized data from various historic pressure computations. The pressure data are 
stored in separate `.csv` files, one per paper, and further classified into groups within files. These files allow 
the user to read the volume, value of pressure, and the error bar of pressure. The data is read from the tables or 
measured using [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/) from figures. If both numerical and graphical 
data exist, copying the numerical data is preferred over digitize data from figures. When the data are obtained from 
figures, each data point is obtained by two measures. One measurement takes place at the top of the error bar (marker), 
and the other at the bottom. Each file contains detailed information about the origin of data. If the data is measured 
from a figure, the figure with measurement markers is provided. The data files could be accessed by either an ascii text 
editor or a spreadsheet. The historic results provided here are from
- [Metropolis 1953](https://doi.org/10.1063/1.1699114)
- [Alder 1962](https://doi.org/10.1103/PhysRev.127.359)
- [Zollweg 1992](https://doi.org/10.1103/PhysRevB.46.11186)
- [Jaster 1999](https://doi.org/10.1103/PhysRevE.59.2594)
- [Jaster 2004](https://doi.org/10.1016/j.physleta.2004.07.055)
- [Mak 2006](https://doi.org/10.1103/PhysRevE.73.065104)
- [Bernard 2011](https://doi.org/10.1103/PhysRevLett.107.155704)
- [Engel 2013](https://doi.org/10.1103/PhysRevE.87.042134)
- [Qi 2014](https://doi.org/10.1039/C4SM00125G)

## Access
The `.csv` files could be accessed directly by either ascii text editors, such as vi, emacs, etc. They could also be 
imported into WYSIWYG editors such as LibreOffice Calc, Google Doc Sheets, Apple Numbers, and Microsoft Office Excel. 
The WYSIWYG editors can also save the imported `.csv` file as `.csv` file. 

`Python/synopsis` contains two Python program that designed to access the `.csv` data files. `read.py` reads in and parses 
the data files and convert the unit of the data, while `plot.py` make plots for the data. The detail and usage of the Python 
program could be found in their docstring.

## Format
The data files are named as `LastNameOfFirstAuthor_YearOfPublication.csv`. The `.csv` file is stored in the directory 
with the same name. If applicable, the screen of the source of the data is stored in the same directory. In each file, 
all the information are store as fields, except for the numerical values. Each field has a field name and a value, and is 
written as FieldName = Value. The spaces before and after the equal sign are mandatory.

The data file begins with a header containing the following fields: 
- Author: the authors of the first publication containing the data
- Year: the year of publication
- Title: title of the paper containing the data
- Journal reference: journal reference of the paper
- Url: url to the paper, sometimes also referred to as DOI
- Data source: = specifying the exact source of digitized data in the paper. This could be either the label of a 
figure or a table, or a piece of descriptive text.
- Extraction method: the method of extracting the data.
- Comment: all the other things worth noting but not covered by the previous six fields. This field is optional.
- Label: the label of the current data file, used by default as the label for plotting by `plot.py`, conventionally as 
"LastNameOfFirstAuthor YearOfPublication"
Each field occupies the first column in each line.

The data within a file is grouped according to the setup of the computation and the representation of data. Each 
group of data begin with an empty line in `.csv` file, namely a line with only two commas when view by an ascii 
text editor. The empty line is followed by two lines containing five tags specifying the state of this group of data. 
The tags are:
- Sampling method: the algorithm implemented to generate the configurations. 
- Pressure calculation method: formula/method implemented to obtain the numerical value of pressure.
- N: number of disks in the calculation. It could be written either as an integer number or a LaTeX formula, for example: 
870, 1024^2.
- Aspect ratio: the aspect ratio of the system during calculation. It could be written either as an integer number or 
in LaTeX format, for example: 1, \sqrt{3} / 2. If the aspect ratio is not specified in the origin of data, it is written 
as unknown.
- Boundary condition: the boundary condition of the run. It can be either non-periodic or periodic.
- X unit: the unit of the x-coordinate, it takes four possible values: "V / V_0", "v / (2 sigma)^2", "packing fraction", 
"reduced density".
- Y unit: the unit of the y-coordinate, it takes two possible values: "beta P V_0 / N", "beta P (2 sigma)^2".
- Type: indicating whether there is an error bar for the data point. It takes two possible values: "error bar" and "marker".

N and Aspect ratio occupy the first two columns in the line following the empty line, and X unit, Y unit, Type occupy the 
line after. The body of data comes after the tags, appearing in the format of multiple lines and three column. The first 
column is the value of the x-coordinate, and the second and the third column are the value and the error bar. In the case 
when the data is given without error bar, the error bar is half the size of the marker. 

The `.png` file stored in the same directory as the `.csv` file provides the figure from which the digitized data is measured. 
They are obtained by taking screenshot of WebPlotDigitizer. 
