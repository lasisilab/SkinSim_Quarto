"""Verify old/corrected BioSkin parameter pairing independently of row IDs."""
import csv
from decimal import Decimal
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def main():
    root = Path(__file__).resolve().parents[1]
    paths = {
        'previous_network': root / 'new_sample_dataset_package/new_sample_skin_props_network_ready.csv',
        'corrected_network': root / 'bioskin_data/corrected_skin_props_network_ready.csv',
        'previous_physical': root / 'new_sample_dataset_package/new_sample_skin_props_raw_physical_units.csv',
        'corrected_physical': root / 'bioskin_data/corrected_skin_props_raw_physical_units.csv',
    }
    old = np.loadtxt(paths['previous_network'], delimiter=',')
    new = np.loadtxt(paths['corrected_network'], delimiter=',')
    permutation = [0, 1, 2, 4, 3, 5]
    aligned = old[:, permutation]
    old_physical = pd.read_csv(paths['previous_physical'])
    new_physical = pd.read_csv(paths['corrected_physical'])
    assert old.shape == new.shape == (4000, 6)
    assert list(old_physical.columns) == list(new_physical.columns)
    assert old_physical.shape == new_physical.shape == (4000, 5)

    # Match independently of row position using the two exactly preserved slots.
    previous_pairs = [tuple(row) for row in aligned[:, 3:5]]
    corrected_pairs = [tuple(row) for row in new[:, 3:5]]
    assert len(set(previous_pairs)) == len(set(corrected_pairs)) == 4000
    pair_lookup = {key: i for i, key in enumerate(previous_pairs)}
    matched_rows = np.array([pair_lookup[key] for key in corrected_pairs])
    assert np.array_equal(matched_rows, np.arange(4000))

    # Also match using the full six-dimensional vector, allowing only tiny drift.
    distances, neighbors = cKDTree(aligned).query(new, k=2)
    assert np.array_equal(neighbors[:, 0], np.arange(4000))
    row_errors = np.max(np.abs(new - aligned), axis=1)
    physical_errors = np.max(np.abs(new_physical.to_numpy() - old_physical.to_numpy()), axis=1)
    assert (row_errors <= 2e-7).all()
    assert (physical_errors <= 1e-7).all()

    # Decimal differences give exact maxima at the CSV's written precision.
    with paths['previous_network'].open(newline='') as handle:
        old_strings = list(csv.reader(handle))
    with paths['corrected_network'].open(newline='') as handle:
        new_strings = list(csv.reader(handle))
    exact_maxima = [str(max(abs(Decimal(n[j]) - Decimal(o[permutation[j]]))
                            for o, n in zip(old_strings, new_strings))) for j in range(6)]
    report = {
        'row_count': 4000,
        'zero_based_previous_column_permutation': permutation,
        'zero_based_swapped_slots': [3, 4],
        'one_based_swapped_columns': [4, 5],
        'row_pairs_verified_by_unique_exact_slot_values': int((matched_rows == np.arange(4000)).sum()),
        'row_pairs_verified_by_nearest_full_parameter_vector': int((neighbors[:, 0] == np.arange(4000)).sum()),
        'rows_exactly_equal_in_all_six_columns_after_swap': int(np.all(new == aligned, axis=1).sum()),
        'rows_matching_with_network_absolute_tolerance_2e_7': int((row_errors <= 2e-7).sum()),
        'max_absolute_difference_per_network_column_exact_decimal': exact_maxima,
        'max_absolute_difference_per_physical_parameter': (new_physical - old_physical).abs().max().to_dict(),
        'rows_matching_physical_parameters_with_absolute_tolerance_1e_7': int((physical_errors <= 1e-7).sum()),
        'max_distance_to_matching_parameter_vector': float(distances[:, 0].max()),
        'min_distance_to_second_closest_parameter_vector': float(distances[:, 1].min()),
        'first_previous_row': old[0].tolist(),
        'first_corrected_row': new[0].tolist(),
        'files': {key: {'path': str(path.relative_to(root)), 'sha256': sha256(path.read_bytes()).hexdigest()}
                  for key, path in paths.items()},
    }
    out = root / 'output/corrected_validation'
    pd.DataFrame({
        'corrected_row_1based': np.arange(1, 4001),
        'matched_previous_row_1based': matched_rows + 1,
        'max_network_difference_after_swap': row_errors,
        'max_physical_parameter_difference': physical_errors,
        'nearest_full_vector_row_1based': neighbors[:, 0] + 1,
        'nearest_full_vector_distance': distances[:, 0],
        'second_closest_full_vector_distance': distances[:, 1],
    }).to_csv(out / 'parameter_pairing_audit.csv', index=False)
    (out / 'parameter_pairing_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
