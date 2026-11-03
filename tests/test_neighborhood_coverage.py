"""Analytic closed-box controls for source-only receptor coverage."""
import numpy as np
import pytest

from aav_benchmark.neighborhood_coverage import closed_box_distances, chain_selection


def key(chain, residue, atom='CA'):
    return (chain, residue, '', 'ALA', atom)


def test_closed_box_inside_faces_edges_and_corners():
    points=np.array([[0,0,0],[1,1,1],[-1,0,0],[4,0,0],[4,5,0],[4,5,13]])
    result=closed_box_distances(points,[0,0,0],[2,2,2])
    assert result.tolist()==pytest.approx([0,0,0,3,5,13])


def test_whole_chain_retained_and_exact_cutoff_is_inclusive():
    keys=[key('A',1),key('A',2),key('B',1),key('C',1)]
    xyz=np.array([[0,0,0],[100,0,0],[8,0,0],[9,0,0]],float)
    six=chain_selection(keys,xyz,[0,0,0],[2,2,2],6)
    eight=chain_selection(keys,xyz,[0,0,0],[2,2,2],8)
    assert six['selected_source_chains']==['A']
    assert six['selected_heavy_atoms']==2  # distant A2 retained with complete chain
    assert six['closest_omitted_atom_to_box_A']==7
    assert eight['selected_source_chains']==['A','B','C']
    assert eight['closest_omitted_atom_to_box_A'] is None


def test_corner_distance_uses_euclidean_not_axiswise_expansion():
    keys=[key('A',1),key('B',1)]
    # B is 6 A beyond each of two box faces: Euclidean sqrt72 > 8 A.
    selection=chain_selection(keys,[[0,0,0],[7,7,0]],[0,0,0],[2,2,2],8)
    assert selection['selected_source_chains']==['A']
    assert selection['closest_omitted_atom_to_box_A']==pytest.approx(np.sqrt(72))


def test_translated_box_preserves_geometry():
    xyz=np.array([[0,0,0],[5,7,9]],float);shift=np.array([31,-22,7])
    assert closed_box_distances(xyz,[0,0,0],[2,4,6]).tolist()==pytest.approx(
        closed_box_distances(xyz+shift,shift,[2,4,6]).tolist())


def test_nearest_omitted_witness_and_chain_minima():
    keys=[key('A',1),key('B',1),key('B',2),key('C',1)]
    selection=chain_selection(keys,[[0,0,0],[11,0,0],[10,0,0],[20,0,0]],[0,0,0],[2,2,2],8)
    assert selection['nearest_omitted_atom']['atom_key']==list(key('B',2))
    assert selection['minimum_distance_per_source_chain_A']=={'A':0,'B':9,'C':19}


@pytest.mark.parametrize('center,size', [([0,0,float('nan')],[2,2,2]),([0,0,0],[2,0,2]),([0,0],[2,2,2])])
def test_invalid_box_rejected(center,size):
    with pytest.raises(ValueError):closed_box_distances([[0,0,0]],center,size)


def test_duplicate_identity_and_invalid_cutoff_rejected():
    with pytest.raises(ValueError,match='Duplicate'):
        chain_selection([key('A',1),key('A',1)],[[0,0,0],[1,0,0]],[0,0,0],[2,2,2],8)
    for cutoff in [0,-1,float('nan')]:
        with pytest.raises(ValueError):chain_selection([key('A',1)],[[0,0,0]],[0,0,0],[2,2,2],cutoff)
