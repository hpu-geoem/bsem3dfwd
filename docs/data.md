# Data and Response Files

## Data File

The data file contains informantion about frequencies, transmitters and receivers. Lines
prefixed by the character `#` are considered as comments. It consists of three
parts.

### Part 1 - Frequencies

This part lists the frequencies.

```text
First line: <# of frequencies>
Following lines list each frequencies:
  <frequency>
  ...
```

### Part 2 - Transmitters

This part lists the transmitters.

```text
One line: <# of transmitters>
Following lines list each transmitters:
  <x> <y> <z> <azimuth> <dip> <current> <dipole length>
  # <x> <y> <z> is the center of transmitter
  # azimuth is the rotation angle of the transmitter from x axis
  # dip is the rotation angle of the transmitter from x-y plane
```

### Part 3 - Receivers

This part lists the receivers.

```text
One line: <# of receivers>
Following lines list each receivers:
  <x> <y> <z>
```

Here is an example:

```text
# frequencies
2 # number of frequencies
 1.00000E+02
 1.00000E+03

# transmitters
1
# X        Y        Y    Azimuth   Dip   Current  Length
0.0       0.0      100.0   0.0     90.0    1.0     100.0

# receivers
1
#     X           Y           Z
0.0000E+00  0.0000E+00  1.0000E-01
```

## Response file

The response file contains a table with a line for each datum. Each line lists
the frequency, transmitter index, coordinates of receivers, and the datum.

```text
<frequency> <transmitter index> <receiver x> <receiver y> <receiver z> <real part of Er> <image part of Er>
```

Here is an example:

```text
1.000000E+02          0            -5.000000E+02     0.000000E+00    1.000000E-01     -6.104814E-08      2.182119E-08
...
```
