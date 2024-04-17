# labtools

Tools useful in Staller Lab for sequence design and analysis. This project is still very much in development and as such is not hosted anywhere official. It is meant for usage by Staller Lab members and collaborators and may change regularly or have bugs at this stage! I plan to rename it and host it somewhere in the future.

#### [See examples](https://massivejords.github.io/tools/docs/_build/html/example.html)

I am still working on implementing test procedures to help with the development process as I make changes. Please let me know if anything looks strange or does not work as intended.

Additionally, please let me know if the trying to install this inconveniences you. I am still trying to understand dependencies and would like to fix any/all annoyances!

## Installation

You can install the latest version using this command (it might be weird as I make adjustments, but probably still better)

```bash 
pip install https://github.com/massivejords/tools/blob/main/dist/labtools-0.1.9-py3-none-any.whl?raw=true
```

Install the older, "stable" VERSION (has some dependency issues that might make it hard to install, but it verifiably does what it is supposed to do)

```bash 
pip install https://github.com/massivejords/tools/blob/main/dist/labtools-0.0.3-py3-none-any.whl?raw=true
```

## Usage

[View the in progress documentation webpage.](https://massivejords.github.io/tools/docs/_build/html/index.html)

[See examples](https://massivejords.github.io/tools/docs/_build/html/example.html)

[Look at functions and classes (with examples)](https://massivejords.github.io/tools/docs/_build/html/autoapi/index.html)

```python
import labtools
from labtools.adtools import sort
```

## License

`labtools` was created by Jordan Stefani. It is licensed under the terms of the MIT license.


## Credits

`labtools` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
