from . import test_all;

from sys import set_int_max_str_digits;

set_int_max_str_digits(100000);

test_all();