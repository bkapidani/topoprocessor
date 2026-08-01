"""
pybind NgOCC module
"""
from __future__ import annotations
import collections.abc
import netgen.libngpy._meshing
import numpy
import numpy.typing
import pyngcore.pyngcore
import typing
__all__: list[str] = ['ApproxParamType', 'ArcOfCircle', 'Axes', 'Axis', 'BSplineCurve', 'BezierCurve', 'BezierSurface', 'Box', 'COMPOUND', 'COMPSOLID', 'Circle', 'Compound', 'Cone', 'ConnectEdgesToWires', 'Cylinder', 'Dir', 'DirectionalInterval', 'EDGE', 'Edge', 'Ellipse', 'Ellipsoid', 'FACE', 'Face', 'From_PyOCC', 'Fuse', 'Geom2d_Curve', 'Geom_Surface', 'Glue', 'HalfSpace', 'ListOfShapes', 'LoadOCCGeometry', 'MakeFillet', 'MakePolygon', 'MakeThickSolid', 'OCCException', 'OCCGeometry', 'Pipe', 'PipeShell', 'Pnt', 'Prism', 'ResetGlobalShapeProperties', 'Revolve', 'SHAPE', 'SHELL', 'SOLID', 'Segment', 'Sew', 'ShapeContinuity', 'Solid', 'Sphere', 'SplineApproximation', 'SplineInterpolation', 'SplineSurfaceApproximation', 'SplineSurfaceInterpolation', 'TestXCAF', 'ThruSections', 'TopAbs_ShapeEnum', 'TopLoc_Location', 'TopoDS_Shape', 'VERTEX', 'Vec', 'Vertex', 'WIRE', 'Wire', 'WorkPlane', 'X', 'Y', 'Z', 'gp_Ax2', 'gp_Ax2d', 'gp_Dir', 'gp_Dir2d', 'gp_GTrsf', 'gp_Mat', 'gp_Pnt', 'gp_Pnt2d', 'gp_Trsf', 'gp_Vec', 'gp_Vec2d', 'occ_version']
class ApproxParamType:
    """
    Wrapper for Approx_ParametrizationType
    
    Members:
    
      Centripetal
    
      ChordLength
    
      IsoParametric
    """
    Centripetal: typing.ClassVar[ApproxParamType]  # value = <ApproxParamType.Centripetal: 1>
    ChordLength: typing.ClassVar[ApproxParamType]  # value = <ApproxParamType.ChordLength: 0>
    IsoParametric: typing.ClassVar[ApproxParamType]  # value = <ApproxParamType.IsoParametric: 2>
    __members__: typing.ClassVar[dict[str, ApproxParamType]]  # value = {'Centripetal': <ApproxParamType.Centripetal: 1>, 'ChordLength': <ApproxParamType.ChordLength: 0>, 'IsoParametric': <ApproxParamType.IsoParametric: 2>}
    @typing.overload
    def __eq__(self, other: ApproxParamType) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: ApproxParamType) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Axes:
    """
    an OCC coordinate system in 3d
    """
    p: gp_Pnt
    @typing.overload
    def __init__(self, p: gp_Pnt = (0, 0, 0), n: gp_Dir = (0, 0, 1), h: gp_Dir = (1, 0, 0)) -> None:
        ...
    @typing.overload
    def __init__(self, axis: Axis) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: gp_Ax2) -> None:
        ...
class Axis:
    """
    an OCC axis in 3d
    """
    def __init__(self, p: gp_Pnt, d: gp_Dir) -> None:
        ...
class Compound(TopoDS_Shape):
    def __init__(self, shapes: collections.abc.Sequence[TopoDS_Shape], separate_layers: bool = False) -> None:
        """
        Create a compound from a list of shapes. If separate_layers is true, assigns layer indices per input shape.
        """
class DirectionalInterval:
    def __and__(self, arg0: DirectionalInterval) -> DirectionalInterval:
        ...
    def __gt__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> DirectionalInterval:
        ...
    def __lt__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> DirectionalInterval:
        ...
    def __str__(self) -> str:
        ...
class Edge(TopoDS_Shape):
    def Extend(self, point: gp_Pnt, continuity: typing.SupportsInt | typing.SupportsIndex = 1, after: bool = True) -> Edge:
        """
        Extend the edge's underlying curve to a target point with G0/G1/G2 continuity.
        """
    def Split(self, *args) -> ...:
        """
        Splits edge at given parameters. Parameters can either be floating values in (0,1), then edge parametrization is used. Or it can be points, then the projection of these points are used for splitting the edge.
        """
    def Tangent(self, s: typing.SupportsFloat | typing.SupportsIndex) -> gp_Vec:
        """
        tangent vector to curve at parameter 's'
        """
    def Value(self, s: typing.SupportsFloat | typing.SupportsIndex) -> gp_Pnt:
        """
        evaluate curve for parameters 's'
        """
    @typing.overload
    def __init__(self, arg0: TopoDS_Shape) -> None:
        """
        Create an edge from a TopoDS_Shape (must be an edge).
        """
    @typing.overload
    def __init__(self, arg0: ..., arg1: TopoDS_Face) -> None:
        """
        Construct an edge from a 2D parametric curve on a face by lifting it to the face's surface.
        """
    @typing.overload
    def __init__(self, arg0: Vertex, arg1: Vertex) -> None:
        """
        Create a straight edge between two vertices.
        """
    @property
    def end(self) -> gp_Pnt:
        """
        end-point of curve
        """
    @property
    def end_tangent(self) -> gp_Vec:
        """
        tangent at end-point
        """
    @property
    def parameter_interval(self) -> tuple[float, float]:
        """
        parameter interval of curve
        """
    @property
    def partition(self) -> pyngcore.pyngcore.Array_D_S | None:
        """
        Optional edge partition parameters for meshing (array of curve parameters).
        """
    @partition.setter
    def partition(self, arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> None:
        ...
    @property
    def start(self) -> gp_Pnt:
        """
        start-point of curve
        """
    @property
    def start_tangent(self) -> gp_Vec:
        """
        tangent at start-point
        """
class Face(TopoDS_Shape):
    def Extend(self, length: typing.SupportsFloat | typing.SupportsIndex, continuity: typing.SupportsInt | typing.SupportsIndex = 1, u_direction: bool = True, after: bool = True) -> Face:
        """
        Extend a bounded face in U or V by a given length with a requested continuity using GeomLib::ExtendSurfByLength. See https://dev.opencascade.org/doc/refman/html/class_geom_lib.html
        """
    def ProjectWire(self, arg0: Wire) -> TopoDS_Shape:
        """
        Project a wire onto a face along the local surface normals using BRepAlgo_NormalProjection. See https://dev.opencascade.org/doc/refman/html/class_b_rep_algo___normal_projection.html
        """
    def WorkPlane(self) -> WorkPlane:
        """
        Create a 2D work plane aligned with the face's surface at (u,v)=(0,0), using the surface normal as the plane normal.
        """
    @typing.overload
    def __init__(self, w: Wire) -> None:
        """
        Create a planar face bounded by a wire.
        """
    @typing.overload
    def __init__(self, f: Face, w: Wire) -> None:
        """
        Create a face on the surface of another face, bounded by a wire.
        """
    @typing.overload
    def __init__(self, f: Face, w: collections.abc.Sequence[Wire]) -> None:
        """
        Create a face on a reference surface and add one or more bounding wires.
        """
    @typing.overload
    def __init__(self, arg0: TopoDS_Shape) -> None:
        """
        Create a face from a TopoDS_Shape (must be a face).
        """
    @property
    def quad_dominated(self) -> bool | None:
        """
        Hint that the face should be meshed with quad-dominated elements.
        """
    @quad_dominated.setter
    def quad_dominated(self, arg1: bool | None) -> None:
        ...
    @property
    def surf(self) -> ...:
        """
        Return the underlying OpenCascade surface of the face.
        """
class Geom2d_Curve:
    def Edge(self) -> Edge:
        """
        Lift the 2D curve to the default plane and return a 3D edge.
        """
    def Face(self) -> Face:
        """
        Create a planar face bounded by the lifted 2D curve.
        """
    def Trim(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> Geom2d_Curve:
        """
        Return a trimmed 2D curve on the parameter interval [u1, u2].
        """
    def Value(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> gp_Pnt2d:
        """
        Evaluate the 2D curve at parameter s.
        """
    def Wire(self) -> Wire:
        """
        Create a wire from the lifted 2D curve on the default plane.
        """
    @property
    def end(self) -> gp_Pnt2d:
        """
        End point of the curve in parameter space.
        """
    @property
    def start(self) -> gp_Pnt2d:
        """
        Start point of the curve in parameter space.
        """
class Geom_Surface:
    def D1(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> tuple[gp_Pnt, gp_Vec, gp_Vec]:
        """
        Evaluate point and first derivatives (du, dv) at parameters (u, v).
        """
    def Normal(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> gp_Dir:
        """
        Compute the surface normal at parameters (u, v) if defined.
        """
    def Value(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> gp_Pnt:
        """
        Evaluate the surface point at parameters (u, v).
        """
class ListOfShapes:
    def Identify(self, other: ListOfShapes, name: str, type: netgen.libngpy._meshing.IdentificationType = ..., trafo: netgen.libngpy._NgOCC.gp_Trsf | netgen.libngpy._NgOCC.gp_GTrsf) -> int:
        """
        Identify shapes for periodic meshing
        """
    def Max(self, dir: gp_Vec) -> typing.Any:
        """
        returns shape where center of gravity is maximal in the direction 'dir'
        """
    def Min(self, dir: gp_Vec) -> typing.Any:
        """
        returns shape where center of gravity is minimal in the direction 'dir'
        """
    @typing.overload
    def Nearest(self, p: gp_Pnt) -> typing.Any:
        """
        returns shape nearest to point 'p'
        """
    @typing.overload
    def Nearest(self, p: gp_Pnt2d) -> typing.Any:
        """
        returns shape nearest to point 'p'
        """
    def Sorted(self, dir: gp_Vec) -> ListOfShapes:
        """
        returns list of shapes, where center of gravity is sorted in direction of 'dir'
        """
    @typing.overload
    def __add__(self, arg0: ListOfShapes) -> ListOfShapes:
        """
        Concatenate two ListOfShapes instances.
        """
    @typing.overload
    def __add__(self, arg0: list) -> ListOfShapes:
        """
        Concatenate a ListOfShapes with a Python list of shapes.
        """
    @typing.overload
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        """
        Return the i-th shape from the list.
        """
    @typing.overload
    def __getitem__(self, arg0: slice) -> ListOfShapes:
        """
        Return a sub-list of shapes using Python slice semantics.
        """
    @typing.overload
    def __getitem__(self, arg0: str) -> ListOfShapes:
        """
        returns list of all shapes named 'name'
        """
    @typing.overload
    def __getitem__(self, arg0: DirectionalInterval) -> ListOfShapes:
        """
        Return shapes whose centers lie inside the given directional interval.
        """
    def __init__(self, arg0: collections.abc.Sequence[TopoDS_Shape]) -> None:
        """
        Create a list of shapes from a Python list.
        """
    def __iter__(self) -> collections.abc.Iterator[typing.Any]:
        """
        Iterate over shapes in the list.
        """
    def __len__(self) -> int:
        """
        Return the number of shapes in the list.
        """
    def __mul__(self, arg0: ListOfShapes) -> ListOfShapes:
        ...
    @property
    def col(self) -> None:
        """
        set col for all elements of list
        """
    @col.setter
    def col(self, arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def edges(self) -> ListOfShapes:
        """
        Return only edge sub-shapes.
        """
    @property
    def faces(self) -> ListOfShapes:
        """
        Return only face sub-shapes.
        """
    @property
    def hpref(self) -> None:
        """
        set hpref for all elements of list
        """
    @hpref.setter
    def hpref(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def maxh(self) -> None:
        """
        set maxh for all elements of list
        """
    @maxh.setter
    def maxh(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def name(self) -> None:
        """
        set name for all elements of list
        """
    @name.setter
    def name(self, arg1: str | None) -> None:
        ...
    @property
    def quad_dominated(self) -> None:
        """
        Set the quad-dominated meshing hint for all shapes in the list.
        """
    @quad_dominated.setter
    def quad_dominated(self, arg1: bool | None) -> None:
        ...
    @property
    def shells(self) -> ListOfShapes:
        """
        Return only shell sub-shapes.
        """
    @property
    def solids(self) -> ListOfShapes:
        """
        Return only solid sub-shapes.
        """
    @property
    def vertices(self) -> ListOfShapes:
        """
        Return only vertex sub-shapes.
        """
    @property
    def wires(self) -> ListOfShapes:
        """
        Return only wire sub-shapes.
        """
class OCCException(Exception):
    pass
class OCCGeometry(netgen.libngpy._meshing.NetgenGeometry):
    """
    Use LoadOCCGeometry to load the geometry from a *.step file.
    """
    def Draw(self) -> None:
        ...
    def GenerateMesh(self, mp: netgen.libngpy._meshing.MeshingParameters = None, comm: pyngcore.pyngcore.MPI_Comm = ..., mesh: netgen.libngpy._meshing.Mesh = None, **kwargs) -> netgen.libngpy._meshing.Mesh:
        """
        Meshing Parameters
        -------------------
        
        maxh: float = 1e10
          Global upper bound for mesh size.
        
        grading: float = 0.3
          Mesh grading how fast the local mesh size can change.
        
        meshsizefilename: str = None
          Load meshsize from file. Can set local mesh size for points
          and along edges. File must have the format:
        
            nr_points
            x1, y1, z1, meshsize
            x2, y2, z2, meshsize
            ...
            xn, yn, zn, meshsize
        
            nr_edges
            x11, y11, z11, x12, y12, z12, meshsize
            ...
            xn1, yn1, zn1, xn2, yn2, zn2, meshsize
        
        segmentsperedge: float = 1.
          Minimal number of segments per edge.
        
        quad_dominated: bool = False
          Quad-dominated surface meshing.
        
        blockfill: bool = True
          Do fast blockfilling.
        
        filldist: float = 0.1
          Block fill up to distance
        
        delaunay: bool = True
          Use delaunay meshing.
        
        delaunay2d : bool = True
          Use delaunay meshing for 2d geometries.
        
        Optimization Parameters
        -----------------------
        
        optimize3d: str = "cmdmustm"
          3d optimization strategy:
            m .. move nodes
            M .. move nodes, cheap functional
            s .. swap faces
            c .. combine elements
            d .. divide elements
            p .. plot, no pause
            P .. plot, Pause
            h .. Histogramm, no pause
            H .. Histogramm, pause
        
        optsteps3d: int = 3
          Number of 3d optimization steps.
        
        optimize2d: str = "smcmSmcmSmcm"
          2d optimization strategy:
            s .. swap, opt 6 lines/node
            S .. swap, optimal elements
            m .. move nodes
            p .. plot, no pause
            P .. plot, pause
            c .. combine
        
        optsteps2d: int = 3
          Number of 2d optimization steps.
        
        elsizeweight: float = 0.2
          Weight of element size w.r.t. element shape in optimization.
        
        
        OCC Specific Meshing Parameters
        -------------------------------
        
        closeedgefac: Optional[float] = 2.
          Factor for meshing close edges, if None it is disabled.
        
        minedgelen: Optional[float] = 0.001
          Minimum edge length to be used for dividing edges to mesh points. If
          None this is disabled.
        """
    def Glue(self) -> None:
        ...
    def Heal(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001, fixsmalledges: bool = True, fixspotstripfaces: bool = True, sewfaces: bool = True, makesolids: bool = True, splitpartitions: bool = False) -> None:
        """
        Heal the OCCGeometry.
        """
    def SetFaceMeshsize(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set maximum meshsize for face fnr. Face numbers are 0 based.
        """
    def __getstate__(self) -> tuple[list]:
        ...
    @typing.overload
    def __init__(self, shape: TopoDS_Shape, dim: typing.SupportsInt | typing.SupportsIndex = 3, copy: bool = False) -> None:
        """
        Create Netgen OCCGeometry from existing TopoDS_Shape
        """
    @typing.overload
    def __init__(self, shape: collections.abc.Sequence[TopoDS_Shape]) -> None:
        """
        Create Netgen OCCGeometry from existing TopoDS_Shape
        """
    @typing.overload
    def __init__(self, filename: str, dim: typing.SupportsInt | typing.SupportsIndex = 3) -> None:
        """
        Load OCC geometry from step, brep or iges file
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def _visualizationData(self) -> dict:
        ...
    @property
    def edges(self) -> ListOfShapes:
        """
        Get edges in order that they will be in the mesh
        """
    @property
    def faces(self) -> ListOfShapes:
        """
        Get faces in order that they will be in the mesh
        """
    @property
    def shape(self) -> TopoDS_Shape:
        ...
    @property
    def solids(self) -> ListOfShapes:
        """
        Get solids in order that they will be in the mesh
        """
    @property
    def vertices(self) -> ListOfShapes:
        """
        Get vertices in order that they will be in the mesh
        """
class ShapeContinuity:
    """
    Wrapper for OCC enum GeomAbs_Shape
    
    Members:
    
      C0
    
      C1
    
      C2
    
      C3
    
      CN
    
      G1
    
      G2
    """
    C0: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.C0: 0>
    C1: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.C1: 2>
    C2: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.C2: 4>
    C3: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.C3: 5>
    CN: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.CN: 6>
    G1: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.G1: 1>
    G2: typing.ClassVar[ShapeContinuity]  # value = <ShapeContinuity.G2: 3>
    __members__: typing.ClassVar[dict[str, ShapeContinuity]]  # value = {'C0': <ShapeContinuity.C0: 0>, 'C1': <ShapeContinuity.C1: 2>, 'C2': <ShapeContinuity.C2: 4>, 'C3': <ShapeContinuity.C3: 5>, 'CN': <ShapeContinuity.CN: 6>, 'G1': <ShapeContinuity.G1: 1>, 'G2': <ShapeContinuity.G2: 3>}
    @typing.overload
    def __eq__(self, other: ShapeContinuity) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: ShapeContinuity) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Solid(TopoDS_Shape):
    def __init__(self, arg0: TopoDS_Shape) -> None:
        """
        Create solid from shell. Shell must consist of topologically closed faces (share vertices and edges).
        """
class TopAbs_ShapeEnum:
    """
    Enumeration of all supported TopoDS_Shapes
    
    Members:
    
      COMPOUND
    
      COMPSOLID
    
      SOLID
    
      SHELL
    
      FACE
    
      WIRE
    
      EDGE
    
      VERTEX
    
      SHAPE
    """
    COMPOUND: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.COMPOUND: 0>
    COMPSOLID: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.COMPSOLID: 1>
    EDGE: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.EDGE: 6>
    FACE: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.FACE: 4>
    SHAPE: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.SHAPE: 8>
    SHELL: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.SHELL: 3>
    SOLID: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.SOLID: 2>
    VERTEX: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.VERTEX: 7>
    WIRE: typing.ClassVar[TopAbs_ShapeEnum]  # value = <TopAbs_ShapeEnum.WIRE: 5>
    __members__: typing.ClassVar[dict[str, TopAbs_ShapeEnum]]  # value = {'COMPOUND': <TopAbs_ShapeEnum.COMPOUND: 0>, 'COMPSOLID': <TopAbs_ShapeEnum.COMPSOLID: 1>, 'SOLID': <TopAbs_ShapeEnum.SOLID: 2>, 'SHELL': <TopAbs_ShapeEnum.SHELL: 3>, 'FACE': <TopAbs_ShapeEnum.FACE: 4>, 'WIRE': <TopAbs_ShapeEnum.WIRE: 5>, 'EDGE': <TopAbs_ShapeEnum.EDGE: 6>, 'VERTEX': <TopAbs_ShapeEnum.VERTEX: 7>, 'SHAPE': <TopAbs_ShapeEnum.SHAPE: 8>}
    @typing.overload
    def __eq__(self, other: TopAbs_ShapeEnum) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: TopAbs_ShapeEnum) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class TopLoc_Location:
    def Transformation(self) -> gp_Trsf:
        ...
    def __init__(self, arg0: gp_Trsf) -> None:
        ...
class TopoDS_Shape:
    def CrossSection(self, plane_axes: Axes) -> TopoDS_Shape:
        """
        Create cross section of shape with plane defined by 'plane_axes' and transfer properties to dim-1 entities
        """
    def Distance(self, arg0: TopoDS_Shape) -> float:
        """
        Compute the minimum distance between two shapes using BRepExtrema_DistShapeShape. See https://dev.opencascade.org/doc/refman/html/class_b_rep_extrema___dist_shape_shape.html
        """
    @typing.overload
    def Extrude(self, h: typing.SupportsFloat | typing.SupportsIndex, dir: netgen.libngpy._NgOCC.gp_Vec | None = None, identify: bool = False, idtype: netgen.libngpy._meshing.IdentificationType = ..., idname: str = 'extrusion') -> TopoDS_Shape:
        """
        extrude shape to thickness 'h', shape must contain a plane surface, optionally give an extrusion direction
        """
    @typing.overload
    def Extrude(self, v: gp_Vec) -> TopoDS_Shape:
        """
        extrude shape by vector 'v'
        """
    def GenerateMesh(self, mp: netgen.libngpy._meshing.MeshingParameters = None, dim: typing.SupportsInt | typing.SupportsIndex = 3, ngs_mesh: bool = True, **kwargs) -> typing.Any:
        """
        Generate a mesh for the shape. Returns an NGSolve mesh if ngs_mesh=True, otherwise a Netgen mesh. Extra keyword arguments are forwarded to OCCGeometry.GenerateMesh.
        """
    def Identify(self, other: TopoDS_Shape, name: str, type: netgen.libngpy._meshing.IdentificationType = ..., trafo: netgen.libngpy._NgOCC.gp_Trsf | netgen.libngpy._NgOCC.gp_GTrsf | None = None) -> int:
        """
        Identify shapes for periodic meshing
        """
    def LimitTolerance(self, tmin: typing.SupportsFloat | typing.SupportsIndex, tmax: typing.SupportsFloat | typing.SupportsIndex = 0.0, type: TopAbs_ShapeEnum = ...) -> None:
        """
        limit tolerance of shape to range [tmin, tmax]
        """
    def Located(self, loc: TopLoc_Location) -> TopoDS_Shape:
        """
        copy shape and sets location of copy
        """
    def MakeChamfer(self, edges: collections.abc.Sequence[TopoDS_Shape], d: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
        """
        make symmetric chamfer for edges 'edges' of distrance 'd'
        """
    @typing.overload
    def MakeFillet(self, fillets: collections.abc.Sequence[tuple[TopoDS_Shape, typing.SupportsFloat | typing.SupportsIndex]]) -> TopoDS_Shape:
        """
        make fillets for shapes of radius 'r'
        """
    @typing.overload
    def MakeFillet(self, edges: collections.abc.Sequence[TopoDS_Shape], r: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
        """
        make fillets for edges 'edges' of radius 'r'
        """
    def MakeThickSolid(self, facestoremove: collections.abc.Sequence[TopoDS_Shape], offset: typing.SupportsFloat | typing.SupportsIndex, tol: typing.SupportsFloat | typing.SupportsIndex, intersection: bool = False, joinType: str = 'arc', removeIntersectingEdges: bool = False) -> TopoDS_Shape:
        """
        makes shell-like solid from faces
        """
    def MakeTriangulation(self) -> None:
        """
        Ensure all faces of the shape have an OpenCascade triangulation (typically via BRepMesh). Useful before querying Poly_Triangulation or exporting to viewers. See https://dev.opencascade.org/doc/refman/html/class_b_rep_mesh___incremental_mesh.html
        """
    @typing.overload
    def Mirror(self, axes: Axes) -> TopoDS_Shape:
        """
        copy shape, and mirror over XY - plane defined by 'axes'
        """
    @typing.overload
    def Mirror(self, axes: Axis) -> TopoDS_Shape:
        """
        copy shape, and rotate by 180 deg around axis 'axis'
        """
    def Move(self, v: gp_Vec) -> typing.Any:
        """
        copy shape, and translate copy by vector 'v'
        """
    def Offset(self, offset: typing.SupportsFloat | typing.SupportsIndex, tol: typing.SupportsFloat | typing.SupportsIndex, intersection: bool = False, joinType: str = 'arc', removeIntersectingEdges: bool = False, identification_name: str | None = None) -> TopoDS_Shape:
        """
        makes shell-like solid from faces
        """
    def Properties(self) -> tuple[typing.Any, typing.Any]:
        """
        returns tuple of shape properties, currently ('mass', 'center'
        """
    def Reversed(self) -> typing.Any:
        """
        Return a copy with the orientation reversed (TopoDS_Shape::Reversed).
        """
    def Revolve(self, axis: Axis, ang: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
        """
        revolve shape around 'axis' by 'ang' degrees
        """
    def Rotate(self, axis: Axis, ang: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
        """
        copy shape, and rotet copy by 'ang' degrees around 'axis'
        """
    def Scale(self, p: gp_Pnt, s: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
        """
        copy shape, and scale copy by factor 's'
        """
    def SetTolerance(self, tol: typing.SupportsFloat | typing.SupportsIndex, stype: TopAbs_ShapeEnum = ...) -> None:
        """
        set (enforce) tolerance of shape to 't'
        """
    def ShapeType(self) -> None:
        """
        deprecated, use 'shape.type' instead
        """
    def SubShapes(self, type: TopAbs_ShapeEnum) -> ...:
        """
        returns list of sub-shapes of type 'type'
        """
    def Triangulation(self) -> ...:
        """
        Extract the face triangulation (Poly_Triangulation) from OpenCascade. If missing, builds it first, then returns the triangle vertex coordinates.
        """
    def UnifySameDomain(self, unifyEdges: bool = True, unifyFaces: bool = True, concatBSplines: bool = True) -> TopoDS_Shape:
        """
        Unify edges and/or faces that lie on the same geometric domain (ShapeUpgrade_UnifySameDomain) and propagate shape properties.
        """
    def WriteBrep(self, filename: str, withTriangles: bool = True, withNormals: bool = False, version: typing.SupportsInt | typing.SupportsIndex | None = None, binary: bool = False) -> None:
        """
        export shape in BREP - format
        """
    def WriteStep(self, filename: str) -> None:
        """
        export shape in STEP - format
        """
    def __add__(self, arg0: TopoDS_Shape) -> TopoDS_Shape:
        """
        fuses shapes
        """
    def __eq__(self, arg0: TopoDS_Shape) -> bool:
        ...
    def __hash__(self) -> int:
        ...
    def __mul__(self, arg0: TopoDS_Shape) -> TopoDS_Shape:
        """
        common of shapes
        """
    def __radd__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> TopoDS_Shape:
        """
        needed for Sum([shapes])
        """
    def __str__(self) -> str:
        ...
    def __sub__(self, arg0: TopoDS_Shape) -> TopoDS_Shape:
        """
        cut of shapes
        """
    def _webgui_data(self) -> dict:
        """
        Return triangulated face/edge data and metadata for web visualization.
        """
    def bc(self, name: str) -> TopoDS_Shape:
        """
        sets 'name' property for all faces of shape
        """
    def mat(self, name: str) -> TopoDS_Shape:
        """
        sets 'name' property to all solids of shape
        """
    @property
    def bounding_box(self) -> tuple[gp_Pnt, gp_Pnt]:
        """
        returns bounding box (pmin, pmax)
        """
    @property
    def center(self) -> gp_Pnt:
        """
        returns center of gravity of shape
        """
    @property
    def col(self) -> typing.Any:
        """
        color of shape as RGB or RGBA - tuple
        """
    @col.setter
    def col(self, arg1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] | None) -> None:
        ...
    @property
    def edges(self) -> ...:
        """
        returns all sub-shapes of type 'EDGE'
        """
    @property
    def faces(self) -> ...:
        """
        returns all sub-shapes of type 'FACE'
        """
    @property
    def hpref(self) -> float:
        """
        number of refinement levels for geometric refinement
        """
    @hpref.setter
    def hpref(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def inertia(self) -> gp_Mat:
        """
        returns matrix of inertia of shape
        """
    @property
    def layer(self) -> int:
        """
        layer of shape
        """
    @layer.setter
    def layer(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def location(self) -> TopLoc_Location:
        """
        Location of shape
        """
    @location.setter
    def location(self, arg1: TopLoc_Location) -> None:
        ...
    @property
    def mass(self) -> float:
        """
        returns mass of shape, what is length, face, or volume
        """
    @property
    def maxh(self) -> float:
        """
        maximal mesh-size for shape
        """
    @maxh.setter
    def maxh(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def name(self) -> str | None:
        """
        'name' of shape
        """
    @name.setter
    def name(self, arg1: str | None) -> None:
        ...
    @property
    def shells(self) -> ...:
        """
        returns all sub-shapes of type 'SHELL'
        """
    @property
    def solids(self) -> ...:
        """
        returns all sub-shapes of type 'SOLID'
        """
    @property
    def type(self) -> TopAbs_ShapeEnum:
        """
        returns type of shape, i.e. 'EDGE', 'FACE', ...
        """
    @property
    def vertices(self) -> ...:
        """
        returns all sub-shapes of type 'VERTEX'
        """
    @property
    def wires(self) -> ...:
        """
        returns all sub-shapes of type 'WIRE'
        """
class Vertex(TopoDS_Shape):
    @typing.overload
    def __init__(self, arg0: TopoDS_Shape) -> None:
        """
        Create a vertex from a TopoDS_Shape (must be a vertex).
        """
    @typing.overload
    def __init__(self, arg0: gp_Pnt) -> None:
        """
        Create a vertex at the given point.
        """
    @property
    def p(self) -> gp_Pnt:
        """
        coordinates of vertex
        """
class Wire(TopoDS_Shape):
    def Offset(self, arg0: TopoDS_Face, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: str, arg3: bool) -> TopoDS_Shape:
        """
        Offset a wire on a supporting face by distance 'dist' with a chosen join type: 'arc' rounds corners with circular arcs, 'tangent' blends with tangent continuity, and 'intersection' keeps sharp corners by intersecting offset segments.
        """
    @typing.overload
    def __init__(self, arg0: Edge) -> None:
        """
        Create a wire from a single edge.
        """
    @typing.overload
    def __init__(self, arg0: collections.abc.Sequence[TopoDS_Shape]) -> None:
        """
        Create a wire from a list of edges and/or wires.
        """
class WorkPlane:
    def Arc(self, r: typing.SupportsFloat | typing.SupportsIndex, ang: typing.SupportsFloat | typing.SupportsIndex, name: str | None = None, maxh: typing.SupportsFloat | typing.SupportsIndex | None = None) -> WorkPlane:
        """
        draw arc tangential to current pos/dir, of radius 'r' and angle 'ang', draw to the left/right if ang is positive/negative
        """
    def ArcTo(self, h: typing.SupportsFloat | typing.SupportsIndex, v: typing.SupportsFloat | typing.SupportsIndex, t: gp_Vec2d, name: str | None = None, maxh: typing.SupportsFloat | typing.SupportsIndex | None = None) -> WorkPlane:
        ...
    @typing.overload
    def Circle(self, h: typing.SupportsFloat | typing.SupportsIndex, v: typing.SupportsFloat | typing.SupportsIndex, r: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        draw circle with center (h,v) and radius 'r'
        """
    @typing.overload
    def Circle(self, r: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        draw circle with center in current position
        """
    def Close(self, name: str | None = None) -> WorkPlane:
        """
        draw line to start point of wire, and finish wire
        """
    def Direction(self, dirh: typing.SupportsFloat | typing.SupportsIndex, dirv: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        reset direction to (dirh, dirv)
        """
    def Ellipse(self, major: typing.SupportsFloat | typing.SupportsIndex, minor: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        draw ellipse with current position as center
        """
    def Face(self) -> Face:
        """
        generate and return face of all wires, resets list of wires
        """
    def Finish(self) -> WorkPlane:
        """
        finish current wire without closing
        """
    def Last(self) -> netgen.libngpy._NgOCC.Wire | None:
        """
        (deprecated) returns current wire
        """
    @typing.overload
    def Line(self, l: typing.SupportsFloat | typing.SupportsIndex, name: str | None = None) -> WorkPlane:
        ...
    @typing.overload
    def Line(self, dx: typing.SupportsFloat | typing.SupportsIndex, dy: typing.SupportsFloat | typing.SupportsIndex, name: str | None = None) -> WorkPlane:
        ...
    def LineTo(self, h: typing.SupportsFloat | typing.SupportsIndex, v: typing.SupportsFloat | typing.SupportsIndex, name: str | None = None) -> WorkPlane:
        """
        draw line to position (h,v)
        """
    def Move(self, l: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        move 'l' from current position and direction, start new wire
        """
    def MoveTo(self, h: typing.SupportsFloat | typing.SupportsIndex, v: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        moveto (h,v), and start new wire
        """
    def NameVertex(self, name: str) -> WorkPlane:
        """
        name vertex at current position
        """
    def Offset(self, d: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        replace current wire by offset curve of distance 'd'
        """
    def Rectangle(self, l: typing.SupportsFloat | typing.SupportsIndex, w: typing.SupportsFloat | typing.SupportsIndex, name: str | None = None) -> WorkPlane:
        """
        draw rectangle, with current position as corner, use current direction
        """
    def RectangleC(self, l: typing.SupportsFloat | typing.SupportsIndex, w: typing.SupportsFloat | typing.SupportsIndex, name: str | None = None) -> WorkPlane:
        """
        draw rectangle, with current position as center, use current direction
        """
    def Reverse(self) -> WorkPlane:
        """
        revert orientation of current wire
        """
    def Rotate(self, ang: typing.SupportsFloat | typing.SupportsIndex) -> WorkPlane:
        """
        rotate current direction by 'ang' degrees
        """
    def Spline(self, points: collections.abc.Sequence[gp_Pnt2d], periodic: bool = False, tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08, tangents: collections.abc.Mapping[typing.SupportsInt | typing.SupportsIndex, gp_Vec2d] = {}, start_from_localpos: bool = True, name: str | None = None) -> WorkPlane:
        """
        draw spline (default: starting from current position, which is implicitly added to given list of points), tangents can be specified for each point (0 refers to starting point)
        """
    def Wire(self) -> netgen.libngpy._NgOCC.Wire | None:
        """
        returns current wire
        """
    def Wires(self) -> ListOfShapes:
        """
        returns all wires
        """
    def __init__(self, axes: Axes = ..., pos: gp_Ax2d = ...) -> None:
        ...
    @property
    def cur_dir(self) -> gp_Vec2d:
        ...
    @property
    def cur_loc(self) -> gp_Pnt2d:
        ...
    @property
    def start_pnt(self) -> gp_Pnt2d:
        ...
class gp_Ax2:
    @typing.overload
    def __init__(self, arg0: gp_Pnt, arg1: gp_Dir) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: gp_Ax3) -> None:
        ...
class gp_Ax2d:
    """
    2d OCC coordinate system
    """
    def __init__(self, p: gp_Pnt2d = (0, 0), d: gp_Dir2d = ...) -> None:
        ...
class gp_Dir:
    """
    3d OCC direction
    """
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> float:
        ...
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    @typing.overload
    def __init__(self, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: gp_Vec) -> None:
        ...
    def __str__(self) -> str:
        ...
class gp_Dir2d:
    """
    2d OCC direction
    """
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    @typing.overload
    def __init__(self, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class gp_GTrsf:
    def __call__(self, arg0: TopoDS_Shape) -> TopoDS_Shape:
        ...
    def __init__(self, mat: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], vec: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = [0.0, 0.0, 0.0]) -> None:
        ...
class gp_Mat:
    """
    3d OCC matrix
    """
    def __getitem__(self, arg0: tuple[typing.SupportsInt | typing.SupportsIndex, typing.SupportsInt | typing.SupportsIndex]) -> float:
        ...
class gp_Pnt:
    """
    3d OCC point
    """
    def __add__(self, arg0: gp_Vec) -> gp_Pnt:
        ...
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> float:
        ...
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    @typing.overload
    def __init__(self, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    @typing.overload
    def __sub__(self, arg0: gp_Pnt) -> gp_Vec:
        ...
    @typing.overload
    def __sub__(self, arg0: gp_Vec) -> gp_Pnt:
        ...
    @property
    def x(self) -> float:
        ...
    @x.setter
    def x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def y(self) -> float:
        ...
    @y.setter
    def y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def z(self) -> float:
        ...
    @z.setter
    def z(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class gp_Pnt2d:
    """
    2d OCC point
    """
    def __add__(self, arg0: gp_Vec2d) -> gp_Pnt2d:
        ...
    @typing.overload
    def __init__(self, arg0: tuple[typing.SupportsFloat | typing.SupportsIndex, typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @typing.overload
    def __init__(self, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    @typing.overload
    def __sub__(self, arg0: gp_Pnt2d) -> gp_Vec2d:
        ...
    @typing.overload
    def __sub__(self, arg0: gp_Vec2d) -> gp_Pnt2d:
        ...
    @property
    def x(self) -> float:
        ...
    @x.setter
    def x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def y(self) -> float:
        ...
    @y.setter
    def y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class gp_Trsf:
    @staticmethod
    def Mirror(arg0: Axis) -> gp_Trsf:
        ...
    @staticmethod
    @typing.overload
    def Rotation(arg0: Axis, arg1: typing.SupportsFloat | typing.SupportsIndex) -> gp_Trsf:
        ...
    @staticmethod
    @typing.overload
    def Rotation(arg0: gp_Pnt, arg1: gp_Dir, arg2: typing.SupportsFloat | typing.SupportsIndex) -> gp_Trsf:
        ...
    @staticmethod
    def Scale(arg0: gp_Pnt, arg1: typing.SupportsFloat | typing.SupportsIndex) -> gp_Trsf:
        ...
    @staticmethod
    @typing.overload
    def Transformation(arg0: Axes) -> gp_Trsf:
        ...
    @staticmethod
    @typing.overload
    def Transformation(arg0: Axes, arg1: Axes) -> gp_Trsf:
        ...
    @staticmethod
    def Translation(arg0: gp_Vec) -> gp_Trsf:
        ...
    def Inverted(self) -> gp_Trsf:
        ...
    def SetMirror(self, arg0: Axis) -> gp_Trsf:
        ...
    def __call__(self, arg0: TopoDS_Shape) -> TopoDS_Shape:
        ...
    def __init__(self) -> None:
        ...
    def __mul__(self, arg0: gp_Trsf) -> gp_Trsf:
        ...
    def __str__(self) -> str:
        ...
class gp_Vec:
    """
    3d OCC vector
    """
    def Norm(self) -> float:
        ...
    def __add__(self, arg0: gp_Vec) -> gp_Vec:
        ...
    def __ge__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> ...:
        ...
    def __gt__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> ...:
        ...
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    @typing.overload
    def __init__(self, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: gp_Dir) -> None:
        ...
    def __le__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> ...:
        ...
    def __lt__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> ...:
        ...
    def __mul__(self, arg0: gp_Vec) -> float:
        ...
    def __neg__(self) -> gp_Vec:
        ...
    def __repr__(self) -> str:
        ...
    def __rmul__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> gp_Vec:
        ...
    def __str__(self) -> str:
        ...
    def __sub__(self, arg0: gp_Vec) -> gp_Vec:
        ...
    def __xor__(self, arg0: gp_Vec) -> gp_Vec:
        ...
    @property
    def x(self) -> float:
        ...
    @x.setter
    def x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def y(self) -> float:
        ...
    @y.setter
    def y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def z(self) -> float:
        ...
    @z.setter
    def z(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class gp_Vec2d:
    """
    2d OCC vector
    """
    def __add__(self, arg0: gp_Vec2d) -> gp_Vec2d:
        ...
    @typing.overload
    def __init__(self, arg0: tuple) -> None:
        ...
    @typing.overload
    def __init__(self, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __neg__(self) -> gp_Vec2d:
        ...
    def __repr__(self: gp_Vec) -> str:
        ...
    def __rmul__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> gp_Vec2d:
        ...
    def __str__(self: gp_Vec) -> str:
        ...
    def __sub__(self, arg0: gp_Vec2d) -> gp_Vec2d:
        ...
    def __xor__(self, arg0: gp_Vec2d) -> float:
        ...
    @property
    def x(self) -> float:
        ...
    @x.setter
    def x(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def y(self) -> float:
        ...
    @y.setter
    def y(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
@typing.overload
def ArcOfCircle(p1: gp_Pnt, p2: gp_Pnt, p3: gp_Pnt) -> Edge:
    """
    Create a circular arc from p1 through p2 to p3.
    """
@typing.overload
def ArcOfCircle(p1: gp_Pnt, v: gp_Vec, p2: gp_Pnt) -> Edge:
    """
    Create a circular arc from p1 to p2 with tangent vector v at p1.
    """
def BSplineCurve(poles: collections.abc.Sequence[gp_Pnt], degree: typing.SupportsInt | typing.SupportsIndex) -> Edge:
    """
    Create a B-spline edge from control points and degree (experimental).
    """
def BezierCurve(points: collections.abc.Sequence[gp_Pnt]) -> Edge:
    """
    Create a Bezier curve from control points.
    """
def BezierSurface(poles: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], weights: typing.Annotated[numpy.typing.ArrayLike, numpy.float64] | None = None, tol: typing.SupportsFloat | typing.SupportsIndex = 1e-07) -> Face:
    """
    Creates a rational Bezier surface with the set of poles and the set of weights. The weights are defaulted to all being 1. If all the weights are identical the surface is considered as non rational. Raises ConstructionError if the number of poles in any direction is greater than MaxDegree + 1 or lower than 2 or CurvePoles and CurveWeights have not the same length or one weight value is lower or equal to Resolution. Returns an occ face with the given tolerance.
    """
def Box(p1: gp_Pnt, p2: gp_Pnt) -> Solid:
    """
    Create a box defined by two opposite corner points p1 and p2.
    """
@typing.overload
def Circle(c: gp_Pnt2d, r: typing.SupportsFloat | typing.SupportsIndex) -> Geom2d_Curve:
    """
    Create a 2D circle curve with center c and radius r.
    """
@typing.overload
def Circle(center: gp_Pnt, normal: gp_Dir, radius: typing.SupportsFloat | typing.SupportsIndex) -> Edge:
    """
    Create a circular edge defined by center, normal, and radius.
    """
def Cone(axis: gp_Ax2, r1: typing.SupportsFloat | typing.SupportsIndex, r2: typing.SupportsFloat | typing.SupportsIndex, h: typing.SupportsFloat | typing.SupportsIndex, angle: typing.SupportsFloat | typing.SupportsIndex) -> Solid:
    """
    Create a cone defined by axis, bottom radius r1, top radius r2, height h, and semi-angle.
    """
def ConnectEdgesToWires(edges: collections.abc.Sequence[TopoDS_Shape], tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08, shared: bool = True) -> list[Wire]:
    """
    Connect edges into one or more wires using ShapeAnalysis_FreeBounds::ConnectEdgesToWires.
    """
@typing.overload
def Cylinder(p: gp_Pnt, d: gp_Dir, r: typing.SupportsFloat | typing.SupportsIndex, h: typing.SupportsFloat | typing.SupportsIndex, bottom: str | None = None, top: str | None = None, mantle: str | None = None) -> typing.Any:
    """
    Create a cylinder at base point p with axis direction d, radius r, and height h. Optional face names: bottom/top/mantle.
    """
@typing.overload
def Cylinder(axis: gp_Ax2, r: typing.SupportsFloat | typing.SupportsIndex, h: typing.SupportsFloat | typing.SupportsIndex) -> Solid:
    """
    Create a cylinder defined by axis, radius, and height.
    """
@typing.overload
def Dir(x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex) -> gp_Dir2d:
    """
    create 2d OCC direction
    """
@typing.overload
def Dir(x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex) -> gp_Dir:
    """
    create 3d OCC direction
    """
@typing.overload
def Dir(d: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> typing.Any:
    """
    create 2d or 3d OCC direction
    """
def Ellipse(axes: gp_Ax2d, major: typing.SupportsFloat | typing.SupportsIndex, minor: typing.SupportsFloat | typing.SupportsIndex) -> Geom2d_Curve:
    """
    Create a 2D ellipse curve defined by axes and major/minor radii.
    """
def Ellipsoid(axes: Axes, r1: typing.SupportsFloat | typing.SupportsIndex, r2: typing.SupportsFloat | typing.SupportsIndex, r3: typing.SupportsFloat | typing.SupportsIndex | None = None) -> TopoDS_Shape:
    """
    Create an ellipsoid aligned with axes, with radii r1, r2, and optional r3 (defaults to r2).
    """
def From_PyOCC(arg0: typing.Any) -> typing.Any:
    """
    Convert a PyOCC SWIG-wrapped TopoDS_Shape into a Netgen TopoDS_Shape view without copying.
    """
def Fuse(arg0: collections.abc.Sequence[TopoDS_Shape]) -> TopoDS_Shape:
    """
    Fuse a list of shapes sequentially (pairwise) using BRepAlgoAPI_Fuse.
    """
@typing.overload
def Glue(shapes: collections.abc.Sequence[TopoDS_Shape]) -> TopoDS_Shape:
    """
    glue together shapes of list
    """
@typing.overload
def Glue(shape: TopoDS_Shape) -> TopoDS_Shape:
    """
    glue together shapes from shape, typically a compound
    """
def HalfSpace(p: gp_Pnt, n: gp_Vec) -> TopoDS_Shape:
    """
    Create a half-space bounded by a plane through point p with normal n.
    """
def LoadOCCGeometry(arg0: os.PathLike | str | bytes) -> netgen.libngpy._meshing.NetgenGeometry:
    ...
def MakeFillet(arg0: TopoDS_Shape, arg1: collections.abc.Sequence[TopoDS_Shape], arg2: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
    """
    deprecated, use 'shape.MakeFillet'
    """
def MakePolygon(verts: collections.abc.Sequence[Vertex]) -> Wire:
    """
    Create a polygonal wire by connecting vertices in order.
    """
def MakeThickSolid(arg0: TopoDS_Shape, arg1: collections.abc.Sequence[TopoDS_Shape], arg2: typing.SupportsFloat | typing.SupportsIndex, arg3: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
    """
    deprecated, use 'shape.MakeThickSolid'
    """
def Pipe(spine: Wire, profile: TopoDS_Shape, twist: tuple[gp_Pnt, typing.SupportsFloat | typing.SupportsIndex] | None = None, auxspine: netgen.libngpy._NgOCC.Wire | None = None) -> TopoDS_Shape:
    """
    Create a pipe by sweeping a profile along a spine wire. If auxspine is provided, uses a pipe shell with the auxiliary spine for orientation.
    """
def PipeShell(spine: Wire, profile: netgen.libngpy._NgOCC.TopoDS_Shape | collections.abc.Sequence[TopoDS_Shape], auxspine: netgen.libngpy._NgOCC.Wire | None = None) -> TopoDS_Shape:
    """
    Create a pipe shell by sweeping one or more profiles along a spine wire. Optionally uses an auxiliary spine to control orientation.
    """
@typing.overload
def Pnt(x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex) -> gp_Pnt2d:
    """
    create 2d OCC point
    """
@typing.overload
def Pnt(x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex) -> gp_Pnt:
    """
    create 3d OCC point
    """
@typing.overload
def Pnt(p: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> typing.Any:
    """
    create 2d or 3d OCC point
    """
def Prism(face: TopoDS_Shape, v: gp_Vec) -> TopoDS_Shape:
    """
    Extrude a face (or shape) along the vector v to create a prism.
    """
def ResetGlobalShapeProperties() -> None:
    """
    Clear cached OpenCascade shape property metadata stored by Netgen.
    """
def Revolve(arg0: TopoDS_Shape, arg1: Axis, arg2: typing.SupportsFloat | typing.SupportsIndex) -> TopoDS_Shape:
    """
    Revolve a shape around an axis by an angle in degrees.
    """
@typing.overload
def Segment(p1: gp_Pnt2d, p2: gp_Pnt2d) -> Geom2d_Curve:
    """
    Create a 2D line segment curve from p1 to p2.
    """
@typing.overload
def Segment(p1: gp_Pnt, p2: gp_Pnt) -> Edge:
    """
    Create a straight edge between two points.
    """
def Sew(faces: collections.abc.Sequence[TopoDS_Shape], tolerance: typing.SupportsFloat | typing.SupportsIndex = 1e-06, non_manifold: bool = False) -> TopoDS_Shape:
    """
    Stitch a list of faces into one or more connected shells.
    
    Parameters
    ----------
    faces : list[TopoDS_Shape]
        Faces or other shapes to sew together.
    tolerance : float, default=1e-6
        Geometric tolerance for merging edges and vertices.
    non_manifold : bool, default=False
        If True, allows edges shared by more than two faces (may produce
        multiple shells). If False, creates only manifold shells suitable
        for solids.
    
    Returns
    -------
    TopoDS_Shape
        The sewed shape containing one or more shells.
    """
def Sphere(c: gp_Pnt, r: typing.SupportsFloat | typing.SupportsIndex) -> Solid:
    """
    Create a sphere with center c and radius r.
    """
@typing.overload
def SplineApproximation(points: collections.abc.Sequence[gp_Pnt2d], approx_type: ApproxParamType = ..., deg_min: typing.SupportsInt | typing.SupportsIndex = 3, deg_max: typing.SupportsInt | typing.SupportsIndex = 8, continuity: ShapeContinuity = ..., tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08) -> Geom2d_Curve:
    """
    Generate a piecewise continuous spline-curve approximating a list of points in 2d.
    
    Parameters
    ----------
    
    points : List|Tuple[gp_Pnt2d]
      List (or tuple) of gp_Pnt.
    
    approx_type : ApproxParamType
      Assumption on location of parameters wrt points.
    
    deg_min : int
      Minimum polynomial degree of splines
    
    deg_max : int
      Maximum polynomial degree of splines
    
    continuity : ShapeContinuity
      Continuity requirement on the approximating surface
    
    tol : float
      Tolerance for the distance from individual points to the approximating curve.
    """
@typing.overload
def SplineApproximation(points: collections.abc.Sequence[gp_Pnt], approx_type: ApproxParamType = ..., deg_min: typing.SupportsInt | typing.SupportsIndex = 3, deg_max: typing.SupportsInt | typing.SupportsIndex = 8, continuity: ShapeContinuity = ..., tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08) -> Edge:
    """
    Generate a piecewise continuous spline-curve approximating a list of points in 3d.
    
    Parameters
    ----------
    
    points : List[gp_Pnt] or Tuple[gp_Pnt]
      List (or tuple) of gp_Pnt.
    
    approx_type : ApproxParamType
      Assumption on location of parameters wrt points.
    
    deg_min : int
      Minimum polynomial degree of splines
    
    deg_max : int
      Maximum polynomial degree of splines
    
    continuity : ShapeContinuity
      Continuity requirement on the approximating surface
    
    tol : float
      Tolerance for the distance from individual points to the approximating curve.
    """
@typing.overload
def SplineInterpolation(points: collections.abc.Sequence[gp_Pnt2d], periodic: bool = False, tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08, tangents: collections.abc.Mapping[typing.SupportsInt | typing.SupportsIndex, gp_Vec2d] = {}) -> Geom2d_Curve:
    """
    Generate a piecewise continuous spline-curve interpolating a list of points in 2d.
    
    Parameters
    ----------
    
    points : List|Tuple[gp_Pnt2d]
      List (or tuple) of gp_Pnt2d.
    
    periodic : bool
      Whether the result should be periodic
    
    tol : float
      Tolerance for the distance between points.
    
    tangents : Dict[int, gp_Vec2d]
      Tangent vectors for the points indicated by the key value (0-based).
    """
@typing.overload
def SplineInterpolation(points: collections.abc.Sequence[gp_Pnt], periodic: bool = False, tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08, tangents: collections.abc.Mapping[typing.SupportsInt | typing.SupportsIndex, gp_Vec] = {}) -> Edge:
    """
    Generate a piecewise continuous spline-curve interpolating a list of points in 3d.
    
    Parameters
    ----------
    
    points : List|Tuple[gp_Pnt]
      List (or tuple) of gp_Pnt
    
    periodic : bool
      Whether the result should be periodic
    
    tol : float
      Tolerance for the distance between points.
    
    tangents : Dict[int, gp_Vec]
      Tangent vectors for the points indicated by the key value (0-based).
    """
def SplineSurfaceApproximation(points: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], approx_type: ApproxParamType = ..., deg_min: typing.SupportsInt | typing.SupportsIndex = 3, deg_max: typing.SupportsInt | typing.SupportsIndex = 8, continuity: ShapeContinuity = ..., tol: typing.SupportsFloat | typing.SupportsIndex = 0.001, periodic: bool = False, degen_tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08) -> Face:
    """
    Generate a piecewise continuous spline-surface approximating an array of points.
    
    Parameters
    ----------
    
    points : np.ndarray
      Array of points coordinates. The first dimension corresponds to the first surface coordinate point
      index, the second dimension to the second surface coordinate point index. The third dimension refers to physical
      coordinates. Such an array can be generated with code like::
    
          px, py = np.meshgrid(*[np.linspace(0, 1, N)]*2)
          points = np.array([[(px[i,j], py[i,j], px[i,j]*py[i,j]**2) for j in range(N)] for i in range(N)])
    
    approx_type : ApproxParamType
      Assumption on location of parameters wrt points.
    
    deg_min : int
      Minimum polynomial degree of splines
    
    deg_max : int
      Maximum polynomial degree of splines
    
    continuity : ShapeContinuity
      Continuity requirement on the approximating surface
    
    tol : float
      Tolerance for the distance from individual points to the approximating surface.
    
    periodic : bool
      Whether the result should be periodic in the first surface parameter
    
    degen_tol : double
      Tolerance for resolution of degenerate edges
    """
def SplineSurfaceInterpolation(points: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], approx_type: ApproxParamType = ..., periodic: bool = False, degen_tol: typing.SupportsFloat | typing.SupportsIndex = 1e-08) -> Face:
    """
    Generate a piecewise continuous spline-surface interpolating an array of points.
    
    Parameters
    ----------
    
    points : np.ndarray
      Array of points coordinates. The first dimension corresponds to the first surface coordinate point
      index, the second dimension to the second surface coordinate point index. The third dimension refers to physical
      coordinates. Such an array can be generated with code like::
    
          px, py = np.meshgrid(*[np.linspace(0, 1, N)]*2)
          points = np.array([[(px[i,j], py[i,j], px[i,j]*py[i,j]**2) for j in range(N)] for i in range(N)])
    
    approx_type : ApproxParamType
      Assumption on location of parameters wrt points.
    
    periodic : bool
      Whether the result should be periodic in the first surface parameter
    
    degen_tol : double
      Tolerance for resolution of degenerate edges
    """
def TestXCAF(shape: TopoDS_Shape = ...) -> None:
    ...
def ThruSections(wires: collections.abc.Sequence[TopoDS_Shape], solid: bool = True) -> TopoDS_Shape:
    """
    Building a loft. This is a shell or solid passing through a set of sections (wires). First and last sections may be vertices. See https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___thru_sections.html#details
    """
@typing.overload
def Vec(x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex) -> gp_Vec2d:
    """
    create 2d OCC point
    """
@typing.overload
def Vec(x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex) -> gp_Vec:
    """
    create 3d OCC point
    """
@typing.overload
def Vec(v: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> typing.Any:
    """
    create 2d or 3d OCC vector
    """
COMPOUND: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.COMPOUND: 0>
COMPSOLID: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.COMPSOLID: 1>
EDGE: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.EDGE: 6>
FACE: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.FACE: 4>
SHAPE: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.SHAPE: 8>
SHELL: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.SHELL: 3>
SOLID: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.SOLID: 2>
VERTEX: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.VERTEX: 7>
WIRE: TopAbs_ShapeEnum  # value = <TopAbs_ShapeEnum.WIRE: 5>
X: gp_Vec  # value = (1, 0, 0)
Y: gp_Vec  # value = (0, 1, 0)
Z: gp_Vec  # value = (0, 0, 1)
occ_version: str = '7.8.1'
