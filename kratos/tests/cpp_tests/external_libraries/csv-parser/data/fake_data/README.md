Data Description
=================

Files
-------
ints.csv
 * Header row with columns "A", "B", ..., "J"
 * 10 columns; each with the first 100 integers
 * "\r\n" separates records
 * **ints.xlsx** has mean and variance calculations
 
ints_newline_sep.csv
 * Same as ints.csv but with "\n" as record separator (not RFC 4180 compliant)

ints_skipline.csv
 * Same as ints.csv, but with an unnecessary line of metadata after the header