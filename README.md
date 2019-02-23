# Convert .res files using PyCORN (adapted to work with pyinstaller)

## Requirments

* matplotlib

## Features

* convert_res.py is a modified PyCORN utility to generate .png plots from .res files
* It does not accept any arguments;
* When run converts every .res file in the same directory as the convert_res.py script.
* Works with PyInstaller allowing convert_res.py to be bundled and easily distrubuted to people not familiar with terminal/console.

## Example usage

```
python convert_res.py
```

## Known issues

On macOS when creating an executable using PyInstaller an error related to tkinter may be issued. Please refer to this [link](https://github.com/pyinstaller/pyinstaller/issues/3753) for a solution.
