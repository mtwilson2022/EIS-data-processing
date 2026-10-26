from eis_processing import smart_dat_to_csv


def test_module_exposes_public_conversion_api():
    assert hasattr(smart_dat_to_csv, "convert_dat_to_csv")
    assert hasattr(smart_dat_to_csv, "combine_trials")
    assert hasattr(smart_dat_to_csv, "main")
