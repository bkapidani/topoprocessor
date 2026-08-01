"""Adapters from supported in-memory mesher objects to the topology core."""

from ._core import Cell, CellKind, Facet, Mesh


def _sequence(value, context):
    try:
        return list(value)
    except TypeError as error:
        raise TypeError(f"{context} must be an iterable") from error


def _integer(value, context):
    for candidate in (value, getattr(value, "nr", None)):
        if candidate is None:
            continue
        try:
            result = int(candidate)
        except (TypeError, ValueError):
            continue
        if result < 0:
            raise ValueError(f"{context} must be non-negative")
        return result
    raise TypeError(f"{context} must be integer-like")


def _point_coordinates(point):
    coordinates = _sequence(getattr(point, "p", point), "Netgen point")
    if len(coordinates) != 3:
        raise ValueError("Netgen points must have three coordinates")
    return tuple(float(coordinate) for coordinate in coordinates)


def _netgen_vertices(element, node_indices, expected, context):
    vertices = _sequence(getattr(element, "vertices", None), context)
    if len(vertices) != expected:
        raise ValueError(f"{context} has {len(vertices)} nodes; expected {expected}")

    result = []
    for vertex in vertices:
        tag = _integer(vertex, f"{context} node")
        try:
            result.append(node_indices[tag])
        except KeyError as error:
            raise ValueError(f"{context} references unknown node {tag}") from error
    return result


def _netgen_mesh(mesh):
    ngmesh = getattr(mesh, "ngmesh", mesh)
    if callable(ngmesh):
        ngmesh = ngmesh()
    for method in ("Points", "Elements3D", "Elements2D"):
        if not callable(getattr(ngmesh, method, None)):
            raise TypeError(
                "expected a Netgen mesh or an NGSolve mesh exposing ngmesh"
            )
    return ngmesh


def from_netgen(mesh):
    """Create a canonical mesh from an in-memory ``netgen.meshing.Mesh``.

    Netgen point identifiers are normalized to zero-based indices. Volume and
    surface ``index`` values become cell and facet labels, respectively.
    """

    ngmesh = _netgen_mesh(mesh)
    raw_points = _sequence(ngmesh.Points(), "Netgen points")
    points = [_point_coordinates(point) for point in raw_points]
    node_indices = {index + 1: index for index in range(len(points))}

    cells = []
    for element in ngmesh.Elements3D():
        vertices = _sequence(getattr(element, "vertices", None), "Netgen cell")
        if len(vertices) == 4:
            kind = CellKind.tetrahedron
        elif len(vertices) == 8:
            kind = CellKind.hexahedron
        else:
            raise ValueError("Netgen cells must be linear tetrahedra or hexahedra")
        nodes = _netgen_vertices(
            element, node_indices, len(vertices), "Netgen cell"
        )
        label = _integer(element.index, "Netgen cell label")
        cells.append(Cell(kind, label, nodes))

    if not cells:
        raise ValueError("Netgen mesh has no volume cells")

    facet_arity = 3 if cells[0].kind == CellKind.tetrahedron else 4
    facets = []
    for element in ngmesh.Elements2D():
        nodes = _netgen_vertices(
            element, node_indices, facet_arity, "Netgen facet"
        )
        facets.append(Facet(_integer(element.index, "Netgen facet label"), nodes))

    return Mesh(points, cells, facets)


def from_ngsolve(mesh):
    """Create a canonical mesh from an ``ngsolve.Mesh``.

    NGSolve is a view over a Netgen mesh, so this adapter unwraps ``mesh.ngmesh``
    and preserves the same label contract as :func:`from_netgen`.
    """

    if not hasattr(mesh, "ngmesh"):
        raise TypeError("expected an NGSolve mesh exposing ngmesh")
    return from_netgen(mesh)


def _gmsh_model(model):
    candidate = getattr(model, "model", model)
    if not callable(getattr(candidate, "getEntities", None)) or not hasattr(
        candidate, "mesh"
    ):
        raise TypeError("expected gmsh.model or the initialized gmsh module")
    return candidate


def _gmsh_label(model, dimension, entity_tag):
    groups = _sequence(
        model.getPhysicalGroupsForEntity(dimension, entity_tag),
        "Gmsh physical groups",
    )
    if len(groups) > 1:
        raise ValueError(
            f"Gmsh entity ({dimension}, {entity_tag}) belongs to multiple "
            "physical groups"
        )
    return 0 if not groups else _integer(groups[0], "Gmsh physical tag")


def _gmsh_nodes(model):
    node_tags, coordinates, _ = model.mesh.getNodes()
    tags = _sequence(node_tags, "Gmsh node tags")
    coordinates = _sequence(coordinates, "Gmsh node coordinates")
    if len(coordinates) != 3 * len(tags):
        raise ValueError("Gmsh returned an inconsistent node coordinate array")

    points = []
    node_indices = {}
    for index, raw_tag in enumerate(tags):
        tag = _integer(raw_tag, "Gmsh node tag")
        if tag == 0:
            raise ValueError("Gmsh node tags must be positive")
        if tag in node_indices:
            raise ValueError(f"Gmsh returned duplicate node tag {tag}")
        node_indices[tag] = index
        offset = 3 * index
        point = tuple(
            float(value) for value in coordinates[offset : offset + 3]
        )
        points.append(point)
    return points, node_indices


def _gmsh_entities(model, dimension, node_indices):
    supported = {2: (3, False), 3: (4, False), 4: (4, True), 5: (8, True)}
    entities = []
    for raw_dimension, raw_entity_tag in model.getEntities(dimension):
        entity_dimension = _integer(raw_dimension, "Gmsh entity dimension")
        entity_tag = _integer(raw_entity_tag, "Gmsh entity tag")
        if entity_dimension != dimension:
            raise ValueError("Gmsh returned an entity with an unexpected dimension")
        label = _gmsh_label(model, dimension, entity_tag)
        element_types, element_tags, element_nodes = model.mesh.getElements(
            dimension, entity_tag
        )
        element_types = _sequence(element_types, "Gmsh element types")
        element_tags = _sequence(element_tags, "Gmsh element tags")
        element_nodes = _sequence(element_nodes, "Gmsh element nodes")
        if not (len(element_types) == len(element_tags) == len(element_nodes)):
            raise ValueError("Gmsh returned inconsistent element arrays")

        for raw_type, tags_for_type, nodes_for_type in zip(
            element_types, element_tags, element_nodes
        ):
            element_type = _integer(raw_type, "Gmsh element type")
            if element_type not in supported:
                kind = "surface" if dimension == 2 else "volume"
                raise ValueError(f"unsupported Gmsh {kind} element type {element_type}")
            arity, is_cell = supported[element_type]
            if is_cell != (dimension == 3):
                raise ValueError(
                    f"Gmsh element type {element_type} has an unexpected dimension"
                )
            tags_for_type = _sequence(tags_for_type, "Gmsh element tags")
            nodes_for_type = _sequence(nodes_for_type, "Gmsh element nodes")
            if len(nodes_for_type) != arity * len(tags_for_type):
                raise ValueError("Gmsh returned an inconsistent element node array")
            for offset in range(0, len(nodes_for_type), arity):
                nodes = []
                for raw_node in nodes_for_type[offset : offset + arity]:
                    node_tag = _integer(raw_node, "Gmsh element node tag")
                    try:
                        nodes.append(node_indices[node_tag])
                    except KeyError as error:
                        raise ValueError(
                            f"Gmsh element references unknown node {node_tag}"
                        ) from error
                if dimension == 3:
                    cell_kind = (
                        CellKind.tetrahedron
                        if element_type == 4
                        else CellKind.hexahedron
                    )
                    entities.append(Cell(cell_kind, label, nodes))
                else:
                    entities.append(Facet(label, nodes))
    return entities


def from_gmsh(model):
    """Create a canonical mesh from the active in-memory ``gmsh.model``.

    Linear triangles/quadrangles and tetrahedra/hexahedra are supported. Each
    model entity must belong to at most one physical group; its physical tag is
    retained as the canonical label, or zero when the entity is ungrouped.
    """

    model = _gmsh_model(model)
    points, node_indices = _gmsh_nodes(model)
    cells = _gmsh_entities(model, 3, node_indices)
    if not cells:
        raise ValueError("Gmsh model has no volume cells")
    facets = _gmsh_entities(model, 2, node_indices)
    return Mesh(points, cells, facets)
