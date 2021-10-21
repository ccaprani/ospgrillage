"""

"""
import pytest
import ospgrillage as og
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Fixtures
@pytest.fixture
def ref_bridge_properties():
    concrete = og.create_material(type="concrete", code="AS5100-2017", grade="50MPa")
    # define sections
    I_beam_section = og.create_section(A=0.896, J=0.133, Iy=0.213, Iz=0.259,
                                       Ay=0.233, Az=0.58)
    slab_section = og.create_section(A=0.04428,
                                     J=2.6e-4, Iy=1.1e-4, Iz=2.42e-4,
                                     Ay=3.69e-1, Az=3.69e-1, unit_width=True)
    exterior_I_beam_section = og.create_section(A=0.044625,
                                                J=2.28e-3, Iy=2.23e-1,
                                                Iz=1.2e-3,
                                                Ay=3.72e-2, Az=3.72e-2)

    # define grillage members
    I_beam = og.create_member(member_name="Intermediate I-beams", section=I_beam_section, material=concrete)
    slab = og.create_member(member_name="concrete slab", section=slab_section, material=concrete)
    exterior_I_beam = og.create_member(member_name="exterior I beams", section=exterior_I_beam_section,
                                       material=concrete)
    return I_beam, slab, exterior_I_beam, concrete

def test_material_command(ref_bridge_properties):
    I_beam, slab, exterior_I_beam, concrete = ref_bridge_properties
    # construct grillage model
    example_bridge = og.OspGrillage(bridge_name="SuperT_10m", long_dim=10, width=7, skew=-42,
                                    num_long_grid=7, num_trans_grid=5, edge_beam_dist=1, mesh_type="Ortho")

    # set grillage member to element groups of grillage model
    example_bridge.set_member(I_beam, member="interior_main_beam")
    example_bridge.set_member(exterior_I_beam, member="exterior_main_beam_1")
    example_bridge.set_member(exterior_I_beam, member="exterior_main_beam_2")
    example_bridge.set_member(exterior_I_beam, member="edge_beam")
    example_bridge.set_member(slab, member="transverse_slab")
    example_bridge.set_member(exterior_I_beam, member="start_edge")
    example_bridge.set_member(exterior_I_beam, member="end_edge")

    example_bridge.create_osp_model(pyfile=False)

    example_bridge.material_command_list