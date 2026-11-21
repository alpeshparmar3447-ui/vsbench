import numpy as np
import pytest

from vsbench.redock import RedockResult, is_success, select_template, worst_engine_rmsd


def _r(name, ad, gn, ok=True):
    return RedockResult(template=name, ok=ok, ad_lowest_rmsd=ad, gn_top_rmsd=gn)


def test_worst_engine_rmsd_is_the_larger_of_the_two():
    assert worst_engine_rmsd(_r("a", 0.8, 1.9)) == 1.9
    assert worst_engine_rmsd(_r("a", 3.0, 1.0)) == 3.0


def test_selection_minimises_the_worst_engine_not_the_average():
    results = [_r("1aaa", 0.5, 4.0), _r("2bbb", 1.4, 1.5), _r("3ccc", 0.1, 3.0)]
    assert select_template(results).template == "2bbb"


def test_selection_ties_break_by_template_name_and_ignore_failures():
    results = [_r("9zzz", 1.0, 1.0), _r("1aaa", 1.0, 1.0), _r("0bad", 0.1, 0.1, ok=False)]
    assert select_template(results).template == "1aaa"


def test_selection_never_uses_nan_results_and_errors_when_nothing_is_valid():
    results = [_r("1aaa", np.nan, 1.0), _r("2bbb", 2.5, 2.6)]
    assert select_template(results).template == "2bbb"
    with pytest.raises(ValueError, match="no template"):
        select_template([_r("x", 1.0, 1.0, ok=False)])


def test_success_threshold_is_strict_and_requires_both_engines():
    assert is_success(_r("a", 1.99, 1.0))
    assert not is_success(_r("a", 2.0, 1.0))
    assert not is_success(_r("a", 0.5, 2.5))
    assert not is_success(_r("a", 0.5, 0.5, ok=False))


def test_selection_uses_no_label_information():
    # The function signature accepts redocking results only; there is nowhere for
    # active/inactive information to enter. Guard that this stays true.
    import inspect
    params = inspect.signature(select_template).parameters
    assert list(params) == ["results"]
