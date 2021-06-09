import pytest
from OpsGrillage import *


# create reference bridge model
@pytest.fixture
def create_reference_bridge():
    # define material
    concrete = UniAxialElasticMaterial(mat_type="Concrete01", mat_vec=[-6.0, -0.004, -6.0, -0.014])

    # define sections
    I_beam_section = Section(op_section_type="Elastic", op_ele_type="elasticBeamColumn", A=0.896, E=3.47E+10,
                             G=2.00E+10,
                             J=0.133, Iy=0.213, Iz=0.259,
                             Ay=0.233, Az=0.58)
    slab_section = Section(op_ele_type="elasticBeamColumn", A=0.04428, E=3.47E+10, G=2.00E+10,
                           J=2.6e-4, Iy=1.1e-4, Iz=2.42e-4,
                           Ay=3.69e-1, Az=3.69e-1, unit_width=True)
    exterior_I_beam_section = Section(op_section_type="Elastic", op_ele_type="elasticBeamColumn", A=0.044625,
                                      E=3.47E+10,
                                      G=2.00E+10, J=2.28e-3, Iy=2.23e-1,
                                      Iz=1.2e-3,
                                      Ay=3.72e-2, Az=3.72e-2)

    # define grillage members
    I_beam = GrillageMember(member_name="Intermediate I-beams", section=I_beam_section, material=concrete)
    slab = GrillageMember(member_name="concrete slab", section=slab_section, material=concrete)
    exterior_I_beam = GrillageMember(member_name="exterior I beams", section=exterior_I_beam_section, material=concrete)

    # construct grillage model
    example_bridge = OpsGrillage(bridge_name="SuperT_10m", long_dim=10, width=7, skew=-42,
                                 num_long_grid=7, num_trans_grid=5, edge_beam_dist=1, mesh_type="Ortho")
    pyfile = False
    example_bridge.create_ops(pyfile=pyfile)

    # set material to grillage
    example_bridge.set_material(concrete)

    # set grillage member to element groups of grillage model
    example_bridge.set_member(I_beam, member="interior_main_beam")
    example_bridge.set_member(exterior_I_beam, member="exterior_main_beam_1")
    example_bridge.set_member(exterior_I_beam, member="exterior_main_beam_2")
    example_bridge.set_member(exterior_I_beam, member="edge_beam")
    example_bridge.set_member(slab, member="transverse_slab")
    example_bridge.set_member(exterior_I_beam, member="start_edge")
    example_bridge.set_member(exterior_I_beam, member="end_edge")

    # Patch load - lane loading
    lane_point_1 = LoadPoint(5, 0, 3, 5)
    lane_point_2 = LoadPoint(8, 0, 3, 5)
    lane_point_3 = LoadPoint(8, 0, 5, 5)
    lane_point_4 = LoadPoint(5, 0, 5, 5)
    Lane = PatchLoading("Lane 1", point1=lane_point_1, point2=lane_point_2, point3=lane_point_3, point4=lane_point_4)

    return example_bridge


# ------------------------------------------------------------------------------------------------------------------
# tests for mesh class object
def test_mesh_negative_skew():
    pass


def test_mesh_positive_skew():
    pass


def test_mesh_skew():
    pass


# ------------------------------------------------------------------------------------------------------------------
# tests for load casees and load setting functions
def test_point_load_getter(create_reference_bridge):  # test get_point_load_nodes() function
    # test if setter and getter is correct              # and assign_point_to_node() function
    ops.wipe()
    example_bridge = create_reference_bridge
    # create reference point load
    location = LoadPoint(5, 0, 2, 20)
    Single = PointLoad(name="single point", point1=location)
    ULS_DL = LoadCase(name="Point")
    ULS_DL.add_load_groups(Single)  # ch
    example_bridge.add_load_case(ULS_DL)
    assert example_bridge.global_load_str == [['ops.load(12, *[0, 1.43019976273292, 0, 0, 0])\n', 'ops.load(17, *[0, '
                                                                                                  '2.56980023726709, '
                                                                                                  '0, 0, 0])\n',
                                               'ops.load(18, *[0, 10.2792009490683, 0, 0, 0])\n', 'ops.load(13, *[0, '
                                                                                                  '5.72079905093166, '
                                                                                                  '0, 0, 0])\n']]


def test_line_load(create_reference_bridge):
    ops.wipe()
    example_bridge = create_reference_bridge
    # create reference line load
    barrierpoint_1 = LoadPoint(3, 0, 3, 2)
    barrierpoint_2 = LoadPoint(10, 0, 3, 2)
    Barrier = LineLoading("Barrier curb load", point1=barrierpoint_1, point2=barrierpoint_2)
    ULS_DL = LoadCase(name="Barrier")
    ULS_DL.add_load_groups(Barrier)  # ch
    example_bridge.add_load_case(ULS_DL)
    assert example_bridge.global_line_int_dict[0] == {7: [[3, 3], [3.1514141550424406, 3.0]],
                                                      8: [[3.1514141550424406, 3.0],
                                                          [4.276919210414739, 2.9999999999999996]],
                                                      11: [[4.276919210414739, 2.9999999999999996],
                                                           [5.402424265787039, 2.9999999999999996]],
                                                      16: [[5.402424265787039, 2.9999999999999996],
                                                           [6.302828310084881, 3.0]], 22: [[6.302828310084881, 3.0],
                                                                                           [7.2271212325636585,
                                                                                            2.9999999999999996]],
                                                      56: [[7.2271212325636585, 2.9999999999999996],
                                                           [8.151414155042438, 2.9999999999999996]],
                                                      62: [[8.151414155042438, 2.9999999999999996],
                                                           [9.075707077521221, 3.0000000000000004]],
                                                      32: [[9.075707077521221, 3.0000000000000004],
                                                           [10.0, 3.0000000000000013]],
                                                      31: [[10.0, 3.0000000000000013], [10, 3]]}


def test_patch_load(create_reference_bridge):
    ops.wipe()
    pass