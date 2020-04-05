import warnings

import numpy as np
import pytest

warnings.filterwarnings("ignore")


@pytest.fixture(scope="session")
def synth():
    from metriclabeling.synthetic import synthetic_dataset

    return synthetic_dataset(n=220, k=12, seed=3)


@pytest.fixture
def rng():
    return np.random.default_rng(0)
