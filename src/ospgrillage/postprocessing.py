# -*- coding: utf-8 -*-
"""
This module contains functions and classes related to post processing processes.
The post processing module is an addition to the currently available post processing
module of OpenSeesPy - this module fills in gaps to
* create envelope from xarray DataSet
* plot force and deflection diagrams from xarray DataSets
"""

import matplotlib.pyplot as plt
import opsvis as opsv
import numpy as np
from typing import TYPE_CHECKING, Union

# if TYPE_CHECKING:
from ospgrillage.load import ShapeFunction
from ospgrillage.utils import solve_zeta_eta

__all__ = [
    "Envelope",
    "PostProcessor",
    "create_envelope",
    "plot_force",
    "plot_defo",
    "plot_bmd",
    "plot_sfd",
    "plot_def",
    "plot_model",
]


# Sentinel for auto-generated titles.  ``title=_AUTO`` means "use the
# default"; ``title=None`` means "no title"; ``title="..."`` is a custom
# override.
_AUTO = object()


# ---------------------------------------------------------------------------
# Lazy import for optional plotly dependency
# ---------------------------------------------------------------------------
def _import_plotly():
    """Lazy-import plotly; raise a clear error if missing."""
    try:
        import plotly.graph_objects as go

        return go
    except ImportError:
        raise ImportError(
            "plotly is required for interactive plots. "
            "Install it with: pip install ospgrillage[gui]"
        )


def create_envelope(**kwargs):
    """
    Create an envelope object for post-processing result envelopes.

    The constructor takes an `xarray` DataSet and kwargs for enveloping options.

    :param ds: Result DataSet from :func:`~ospgrillage.osp_grillage.OspGrillage.get_results`.
    :type ds: xarray.Dataset
    :param load_effect: Specific load effect component to envelope.
    :type load_effect: str, optional
    :param array: Data array to envelope — either ``"displacements"`` or ``"forces"``.
    :type array: str, optional
    :param value_mode: If ``True``, return raw envelope values. Defaults to ``True``.
    :type value_mode: bool, optional
    :param query_mode: If ``True``, return the load case coordinates at the envelope extrema. Defaults to ``False``.
    :type query_mode: bool, optional
    :param extrema: Envelope direction — either ``"min"`` or ``"max"``.
    :type extrema: str, optional
    :param elements: Specific element tags to include in the envelope.
    :type elements: list, optional
    :param nodes: Specific node tags to include in the envelope.
    :type nodes: list, optional
    :returns: :class:`Envelope` object.
    """
    return Envelope(**kwargs)


class Envelope:
    """
    Class for defining envelope of :class:`~ospgrillage.osp_grillage.OspGrillage`'s
    `xarray` of result.

    A :func:`Envelope.get` method is provided that returns an enveloped `xarray` based
    on the specified input parameters of this class.


    """

    def __init__(self, ds, load_effect: str = None, **kwargs):
        """

        :param ds: Data set from
                   :func:`~ospgrillage.osp_grillage.OspGrillage.get_results` . note Combination
        :type ds: Xarray
        :param load_effect: Specific load effect to envelope.
        :type load_effect: str
        :param array: Data array to envelope — either ``"displacements"`` or ``"forces"``.
        :type array: str, optional
        :param value_mode: If ``True``, return raw envelope values. Defaults to ``True``.
        :type value_mode: bool, optional
        :param query_mode: If ``True``, return the load case coordinate at the envelope extrema. Defaults to ``False``.
        :type query_mode: bool, optional
        :param extrema: Envelope direction — either ``"min"`` or ``"max"``.
        :type extrema: str, optional
        :param elements: Specific element tags to include in the envelope.
        :type elements: list, optional
        :param nodes: Specific node tags to include in the envelope.
        :type nodes: list, optional

        """
        self.value = True
        self.ds = ds
        if ds is None:
            return

        # instantiate variables
        self.load_effect = (
            load_effect  # array load effect either displacements or forces
        )
        self.envelope_ds = None
        # default xarray function name
        self.xarray_command = {
            "query": ["idxmax", "idxmin"],
            "minmax value": ["max", "min"],
            "index": ["argmax", "argmin"],
        }
        self.selected_xarray_command = []
        # get keyword args
        self.elements = kwargs.get(
            "elements", None
        )  # specific elements to query/envelope
        self.nodes = kwargs.get("nodes", None)  # specific nodes to query/envelope
        self.component = kwargs.get(
            "load_effect", None
        )  # specific load effect to query
        self.array = kwargs.get("array", "displacements")
        self.value_mode = kwargs.get("value_mode", True)
        self.query_mode = kwargs.get("query_mode", False)  # default query mode
        self.extrema = kwargs.get("extrema", "max")

        # check variables
        if self.load_effect is None:
            raise ValueError(
                "Missing argument for load_effect=: Hint requires a namestring of load"
                "effect type based on the Component dimension of the ospgrillage data"
                "set result format"
            )

        # process variables
        self.extrema_index = 0 if self.extrema == "max" else 1  # minima
        if self.query_mode:
            self.selected_xarray_command = self.xarray_command["query"][
                self.extrema_index
            ]
        elif self.value_mode:
            self.selected_xarray_command = self.xarray_command["minmax value"][
                self.extrema_index
            ]
        else:  # default to argmax/ argmin
            self.selected_xarray_command = self.xarray_command["index"][
                self.extrema_index
            ]

        # NOTE: element/component filtering via self.elements and self.component
        # is not yet implemented — get() currently returns the full reduced array.
        # A future enhancement should build a sel() dict from these kwargs and
        # call .sel(**sel_kwargs) in get().

    def get(self):
        """

        :returns: The enveloped xarray DataArray with envelope values applied.
        :rtype: xarray.DataArray
        """
        da = getattr(self.ds, self.array)
        return getattr(da, self.selected_xarray_command)(dim="Loadcase")


# ---------------------------------------------------------------------------
# Force component lookup tables
# ---------------------------------------------------------------------------
_FORCE_COMP_DICT = {"Fx": 0, "Fy": 1, "Fz": 2, "Mx": 3, "My": 4, "Mz": 5}
_FORCE_COMP_FACTOR = {"Fx": 1, "Fy": 1, "Fz": 1, "Mx": 1, "My": 1, "Mz": -1}

_FORCE_COMPONENTS = [
    "Vx_i", "Vy_i", "Vz_i", "Mx_i", "My_i", "Mz_i",
    "Vx_j", "Vy_j", "Vz_j", "Mx_j", "My_j", "Mz_j",
]


# ---------------------------------------------------------------------------
# Data extraction helpers (backend-agnostic)
# ---------------------------------------------------------------------------
def _extract_force_data(ospgrillage_obj, result_obj, component, member,
                        loadcase=None, option="elements"):
    """Extract force distribution data for a single member.

    Iterates over every element in the requested *member* group, calls
    :func:`opsvis.section_force_distribution_3d` for each, and collects the
    results into three parallel lists.

    :param ospgrillage_obj: Grillage model object.
    :type ospgrillage_obj: OspGrillage
    :param result_obj: xarray DataSet of analysis results.
    :param component: Force component name (``"Fx"``, ``"Fy"``, ``"Fz"``,
        ``"Mx"``, ``"My"``, ``"Mz"``) or integer index.
    :param member: Member group name string.
    :param loadcase: Load case name.  ``None`` uses the first load case.
    :param option: Element query option passed to
        :meth:`~ospgrillage.osp_grillage.OspGrillage.get_element`.
    :returns: ``(all_xx, all_zz, all_values)`` — each a list of
        per-element :class:`numpy.ndarray`.  ``all_values`` entries are
        already multiplied by the component sign factor.
    :rtype: tuple[list, list, list]
    """
    comp_index = component if isinstance(component, int) else _FORCE_COMP_DICT[component]
    factor = _FORCE_COMP_FACTOR[component] if isinstance(component, str) else 1

    nodes = ospgrillage_obj.get_nodes()
    eletag = ospgrillage_obj.get_element(member=member, options=option)

    if ospgrillage_obj.model_type == "shell_beam":
        force_result = result_obj.forces_beam
    else:
        force_result = result_obj.forces

    all_xx, all_zz, all_values = [], [], []

    for ele in eletag:
        # get force components for this element
        if loadcase:
            ele_components = force_result.sel(
                Loadcase=loadcase, Element=ele, Component=_FORCE_COMPONENTS,
            ).values
        else:
            ele_components = force_result.sel(
                Element=ele, Component=_FORCE_COMPONENTS,
            )[0].values

        # get nodes of element
        if ospgrillage_obj.model_type == "shell_beam":
            ele_node = result_obj.ele_nodes_shell.sel(Element=ele)
            if all(np.isnan(result_obj.ele_nodes_shell.sel(Element=ele))):
                ele_node = result_obj.ele_nodes_beam.sel(Element=ele)
            ele_node = ele_node[~np.isnan(ele_node)]
        else:
            ele_node = result_obj.ele_nodes.sel(Element=ele)

        xx = [nodes[n]["coordinate"][0] for n in ele_node.values]
        yy = [nodes[n]["coordinate"][1] for n in ele_node.values]
        zz = [nodes[n]["coordinate"][2] for n in ele_node.values]

        # opsvis section force distribution
        s, *_ = opsv.section_force_distribution_3d(
            ecrd=np.column_stack([xx, yy, zz]),
            pl=ele_components,
        )

        all_xx.append(np.array(xx))
        all_zz.append(np.array(zz))
        all_values.append(s[:, comp_index] * factor)

    return all_xx, all_zz, all_values


def _extract_defo_data(ospgrillage_obj, result_obj, member, component="y",
                       loadcase=None, option="nodes"):
    """Extract displacement data for a single member.

    Iterates over the nodes belonging to *member* and collects their
    coordinates and displacement values.

    :param ospgrillage_obj: Grillage model object.
    :type ospgrillage_obj: OspGrillage
    :param result_obj: xarray DataSet of analysis results.
    :param member: Member group name string.
    :param component: Displacement component (``"x"``, ``"y"``, ``"z"``).
        Defaults to ``"y"`` (vertical deflection).
    :param loadcase: Load case name.  ``None`` uses the first load case.
    :param option: Node query option passed to
        :meth:`~ospgrillage.osp_grillage.OspGrillage.get_element`.
    :returns: ``(xx_list, zz_list, values_list)`` — node x-coordinates,
        z-coordinates, and displacement scalars.
    :rtype: tuple[list, list, list]
    """
    plot_option = option if option is not None else "nodes"
    dis_comp = component if component is not None else "y"

    nodes = ospgrillage_obj.get_nodes()
    nodes_to_plot = ospgrillage_obj.get_element(member=member, options=plot_option)[0]

    xx_list, zz_list, values_list = [], [], []

    for node in nodes_to_plot:
        if loadcase:
            disp = result_obj.displacements.sel(
                Component=dis_comp, Node=node, Loadcase=loadcase,
            ).values
        else:
            disp = result_obj.displacements.sel(
                Component=dis_comp, Node=node,
            )[0].values
        xx_list.append(nodes[node]["coordinate"][0])
        zz_list.append(nodes[node]["coordinate"][2])
        values_list.append(float(disp))

    return xx_list, zz_list, values_list


# ---------------------------------------------------------------------------
# Plotly 3D figure builders
# ---------------------------------------------------------------------------
def _plotly_3d_force(ospgrillage_obj, result_obj, component, members,
                     loadcase=None, comp_label=None, *,
                     fig=None, figsize=None, scale=1.0, title=_AUTO,
                     alpha=1.0):
    """Build an interactive 3D Plotly figure of force diagrams.

    Each member is drawn as a ``Scatter3d`` trace with:

    * x = longitudinal position (m)
    * y = transverse position (z-coordinate from the grillage model, m)
    * z = force component value

    Vertical drop-lines connect the diagram to a dotted zero-baseline so
    the diagram shape is visible from any rotation angle.

    :param ospgrillage_obj: Grillage model object.
    :param result_obj: xarray DataSet of analysis results.
    :param component: Force component name (e.g. ``"Mz"``, ``"Fy"``).
    :param members: List of member group name strings to include.
    :param loadcase: Load case name.  ``None`` uses the first load case.
    :param comp_label: Axis / title label.  Defaults to *component*.
    :param fig: Existing Plotly ``Figure`` to add traces to.
    :param figsize: Figure size in inches ``(width, height)``.
    :param scale: Multiply plotted values by this factor.
    :param title: Plot title.  Default auto-generates from component.
        ``None`` suppresses the title.
    :param alpha: Trace opacity (0–1).
    :returns: Interactive 3D figure.
    :rtype: :class:`plotly.graph_objects.Figure`
    """
    go = _import_plotly()
    if fig is None:
        fig = go.Figure()
    label = comp_label or component

    colours = ["black", "blue", "red", "green", "orange"]

    for idx, member in enumerate(members):
        all_xx, all_zz, all_vals = _extract_force_data(
            ospgrillage_obj, result_obj, component, member, loadcase,
        )
        colour = colours[idx % len(colours)]

        for ex, ez, ev in zip(all_xx, all_zz, all_vals):
            sv = ev * scale
            # Force diagram line
            fig.add_trace(go.Scatter3d(
                x=ex, y=ez, z=sv,
                mode="lines",
                line=dict(color=colour, width=4),
                opacity=alpha,
                name=member,
                showlegend=False,
            ))
            # Baseline (zero) line
            fig.add_trace(go.Scatter3d(
                x=ex, y=ez, z=np.zeros_like(sv),
                mode="lines",
                line=dict(color=colour, width=1, dash="dot"),
                opacity=alpha,
                showlegend=False,
            ))
            # Vertical drop-lines connecting diagram to baseline
            for xi, zi, vi in zip(ex, ez, sv):
                fig.add_trace(go.Scatter3d(
                    x=[xi, xi], y=[zi, zi], z=[0, vi],
                    mode="lines",
                    line=dict(color=colour, width=1),
                    opacity=alpha,
                    showlegend=False,
                ))

        # One legend entry per member
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None],
            mode="lines",
            line=dict(color=colour, width=4),
            name=member,
        ))

    layout_kw = dict(
        scene=dict(
            xaxis_title="x (m)",
            yaxis_title="z (m)",
            zaxis_title=label,
        ),
        legend_title="Member",
    )
    if title is _AUTO:
        layout_kw["title"] = f"{label} Diagram"
    elif title is not None:
        layout_kw["title"] = title
    if figsize is not None:
        layout_kw["width"] = figsize[0] * 100
        layout_kw["height"] = figsize[1] * 100

    fig.update_layout(**layout_kw)
    return fig


def _plotly_3d_defo(ospgrillage_obj, result_obj, members, component="y",
                    loadcase=None, *,
                    fig=None, figsize=None, scale=1.0, title=_AUTO):
    """Build an interactive 3D Plotly figure of deflected shapes.

    Each member is drawn as a ``Scatter3d`` trace with markers, accompanied
    by a dotted zero-baseline.

    :param ospgrillage_obj: Grillage model object.
    :param result_obj: xarray DataSet of analysis results.
    :param members: List of member group name strings to include.
    :param component: Displacement component (default ``"y"``).
    :param loadcase: Load case name.  ``None`` uses the first load case.
    :param fig: Existing Plotly ``Figure`` to add traces to.
    :param figsize: Figure size in inches ``(width, height)``.
    :param scale: Multiply plotted values by this factor.
    :param title: Plot title.  Default auto-generates from component.
        ``None`` suppresses the title.
    :returns: Interactive 3D figure.
    :rtype: :class:`plotly.graph_objects.Figure`
    """
    go = _import_plotly()
    if fig is None:
        fig = go.Figure()

    colours = ["blue", "red", "green", "orange", "black"]

    for idx, member in enumerate(members):
        xx, zz, vals = _extract_defo_data(
            ospgrillage_obj, result_obj, member, component, loadcase,
        )
        colour = colours[idx % len(colours)]
        sv = np.array(vals) * scale

        fig.add_trace(go.Scatter3d(
            x=xx, y=zz, z=sv.tolist(),
            mode="lines+markers",
            line=dict(color=colour, width=4),
            marker=dict(size=3, color=colour),
            name=member,
        ))
        # Baseline
        fig.add_trace(go.Scatter3d(
            x=xx, y=zz, z=[0.0] * len(vals),
            mode="lines",
            line=dict(color=colour, width=1, dash="dot"),
            showlegend=False,
        ))

    layout_kw = dict(
        scene=dict(
            xaxis_title="x (m)",
            yaxis_title="z (m)",
            zaxis_title=f"displacement ({component})",
        ),
        legend_title="Member",
    )
    if title is _AUTO:
        layout_kw["title"] = f"Deflection ({component})"
    elif title is not None:
        layout_kw["title"] = title
    if figsize is not None:
        layout_kw["width"] = figsize[0] * 100
        layout_kw["height"] = figsize[1] * 100

    fig.update_layout(**layout_kw)
    return fig


# ---------------------------------------------------------------------------
# Matplotlib plotting functions
# ---------------------------------------------------------------------------
def plot_force(
    ospgrillage_obj,
    result_obj=None,
    component=None,
    member: str = None,
    option: str = "elements",
    loadcase: str = None,
    *,
    figsize=None,
    ax=None,
    scale: float = 1.0,
    title=_AUTO,
    color: str = "k",
    fill: bool = True,
    alpha: float = 0.4,
    show: bool = False,
):
    """
    Plot a force diagram for a grillage model result for a specified component and load case.

    .. note::
        For "shell_beam" model type, the function only plots the force diagrams for beam elements only.

    :param ospgrillage_obj: Grillage model object.
    :type ospgrillage_obj: OspGrillage
    :param result_obj: xarray DataSet of results.
    :type result_obj: xarray DataSet
    :param component: Force component to plot (e.g. ``"Mz"``, ``"Fy"``).
    :type component: str
    :param member: Member group name (required).
    :type member: str
    :param option: Element query option.
    :type option: str
    :param loadcase: Load case name.  If ``None``, uses the first load case.
    :type loadcase: str
    :param figsize: Figure size in inches ``(width, height)``.
        Ignored when *ax* is provided.
    :type figsize: tuple, optional
    :param ax: Existing matplotlib Axes to plot on.  When provided the
        function draws on this axes instead of creating a new figure.
    :type ax: :class:`~matplotlib.axes.Axes`, optional
    :param scale: Multiply plotted values by this factor (e.g. ``0.001``
        to convert N to kN).
    :type scale: float
    :param title: Plot title.  Default (``_AUTO``) uses the member name.
        Pass a string to override, or ``None`` to suppress the title.
    :type title: str or None
    :param color: Line / fill colour.
    :type color: str
    :param fill: If ``True`` (default), shade the area under the diagram.
    :type fill: bool
    :param alpha: Fill transparency (0 = transparent, 1 = opaque).
    :type alpha: float
    :param show: If ``True``, call ``plt.show()`` before returning.
    :type show: bool
    :returns: Matplotlib axes (use ``ax.get_figure()`` to obtain the figure).
    :rtype: :class:`~matplotlib.axes.Axes`
    """
    if member is None:
        raise ValueError("Missing argument: member= is required")

    all_xx, _all_zz, all_values = _extract_force_data(
        ospgrillage_obj, result_obj, component, member, loadcase, option,
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    for xx, vals in zip(all_xx, all_values):
        scaled = vals * scale
        ax.plot(xx, scaled, "-", color=color)
        if fill:
            ax.fill_between(xx, scaled, np.zeros_like(scaled),
                            color=color, alpha=alpha)

    if title is _AUTO:
        ax.set_title(member)
    elif title is not None:
        ax.set_title(title)
    ax.set_xlabel("x (m) ")
    ax.set_ylabel(component)
    fig.tight_layout()

    if show:
        plt.show()

    return ax


def plot_defo(
    ospgrillage_obj,
    result_obj=None,
    member: str = None,
    component: str = None,
    option: str = "nodes",
    loadcase: str = None,
    *,
    figsize=None,
    ax=None,
    scale: float = 1.0,
    title=_AUTO,
    color: str = "b",
    show: bool = False,
):
    """
    Plot a displacement diagram for a grillage model result for a specified component and load case.

    .. note::
        For "shell_beam" model type, the function only plots the force diagrams for beam elements only.

    :param ospgrillage_obj: Grillage model object.
    :type ospgrillage_obj: OspGrillage
    :param result_obj: xarray DataSet of results.
    :type result_obj: xarray DataSet
    :param component: Displacement component (default ``"y"``).
    :type component: str
    :param member: Member group name (required).
    :type member: str
    :param option: Node query option, either ``"nodes"`` or ``"element"``.
    :type option: str
    :param loadcase: Load case name.  If ``None``, uses the first load case.
    :type loadcase: str
    :param figsize: Figure size in inches ``(width, height)``.
        Ignored when *ax* is provided.
    :type figsize: tuple, optional
    :param ax: Existing matplotlib Axes to plot on.  When provided the
        function draws on this axes instead of creating a new figure.
    :type ax: :class:`~matplotlib.axes.Axes`, optional
    :param scale: Multiply plotted values by this factor.
    :type scale: float
    :param title: Plot title.  Default (``_AUTO``) uses the member name.
        Pass a string to override, or ``None`` to suppress the title.
    :type title: str or None
    :param color: Line colour.
    :type color: str
    :param show: If ``True``, call ``plt.show()`` before returning.
    :type show: bool
    :returns: Matplotlib axes (use ``ax.get_figure()`` to obtain the figure).
    :rtype: :class:`~matplotlib.axes.Axes`
    """
    if member is None:
        raise ValueError("Missing argument: member= is required")

    dis_comp = component if component is not None else "y"
    xx_list, _zz_list, values_list = _extract_defo_data(
        ospgrillage_obj, result_obj, member, dis_comp, loadcase, option,
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    scaled = np.array(values_list) * scale
    ax.plot(xx_list, scaled, "-", color=color)

    if title is _AUTO:
        ax.set_title(member)
    elif title is not None:
        ax.set_title(title)
    ax.set_xlabel("x (m) ")
    ax.set_ylabel(dis_comp)
    fig.tight_layout()

    if show:
        plt.show()

    return ax


# ---------------------------------------------------------------------------
# Convenience plotting wrappers
# ---------------------------------------------------------------------------
_MAIN_BEAM_MEMBERS = [
    "exterior_main_beam_1",
    "interior_main_beam",
    "exterior_main_beam_2",
]


def plot_bmd(ospgrillage_obj, result_obj=None, member=None, loadcase=None,
             backend="matplotlib", **kwargs):
    """
    Plot bending moment diagram (Mz) for one or all main beams.

    When *member* is ``None``, iterates over the three main-beam member groups
    and returns a list of figures.

    :param ospgrillage_obj: Grillage model object.
    :param result_obj: xarray DataSet of results.
    :param member: Member name. If ``None``, plots all main beams.
    :param loadcase: Load case name. If ``None``, uses the first load case.
    :param backend: ``"matplotlib"`` (default, static) or ``"plotly"`` (interactive 3D).
    :param \\**kwargs: Forwarded to the underlying renderer.  See
        :func:`plot_force` (matplotlib) or the Plotly builder for accepted
        keyword arguments such as *figsize*, *ax*, *scale*, *title*,
        *color*, *fill*, *alpha*, and *show*.
    :returns: Single axes when *member* is given, else list of axes.
        For ``backend="plotly"``, returns a single :class:`plotly.graph_objects.Figure`.
    """
    members = [member] if member else _MAIN_BEAM_MEMBERS
    if backend == "plotly":
        plotly_kw = {k: v for k, v in kwargs.items()
                     if k in ("figsize", "scale", "title", "alpha")}
        plotly_kw["fig"] = kwargs.get("ax", None)
        return _plotly_3d_force(
            ospgrillage_obj, result_obj, component="Mz",
            members=members, loadcase=loadcase, comp_label="Mz",
            **plotly_kw,
        )
    # matplotlib path
    if member is not None:
        return plot_force(
            ospgrillage_obj, result_obj, component="Mz",
            member=member, loadcase=loadcase, **kwargs,
        )
    figs = []
    for m in _MAIN_BEAM_MEMBERS:
        try:
            figs.append(
                plot_force(
                    ospgrillage_obj, result_obj, component="Mz",
                    member=m, loadcase=loadcase, **kwargs,
                )
            )
        except (ValueError, KeyError, IndexError):
            pass
    return figs


def plot_sfd(ospgrillage_obj, result_obj=None, member=None, loadcase=None,
             backend="matplotlib", **kwargs):
    """
    Plot shear force diagram (Fy) for one or all main beams.

    When *member* is ``None``, iterates over the three main-beam member groups
    and returns a list of figures.

    :param ospgrillage_obj: Grillage model object.
    :param result_obj: xarray DataSet of results.
    :param member: Member name. If ``None``, plots all main beams.
    :param loadcase: Load case name. If ``None``, uses the first load case.
    :param backend: ``"matplotlib"`` (default, static) or ``"plotly"`` (interactive 3D).
    :param \\**kwargs: Forwarded to the underlying renderer.  See
        :func:`plot_force` for accepted keyword arguments.
    :returns: Single axes when *member* is given, else list of axes.
        For ``backend="plotly"``, returns a single :class:`plotly.graph_objects.Figure`.
    """
    members = [member] if member else _MAIN_BEAM_MEMBERS
    if backend == "plotly":
        plotly_kw = {k: v for k, v in kwargs.items()
                     if k in ("figsize", "scale", "title", "alpha")}
        plotly_kw["fig"] = kwargs.get("ax", None)
        return _plotly_3d_force(
            ospgrillage_obj, result_obj, component="Fy",
            members=members, loadcase=loadcase, comp_label="Fy",
            **plotly_kw,
        )
    # matplotlib path
    if member is not None:
        return plot_force(
            ospgrillage_obj, result_obj, component="Fy",
            member=member, loadcase=loadcase, **kwargs,
        )
    figs = []
    for m in _MAIN_BEAM_MEMBERS:
        try:
            figs.append(
                plot_force(
                    ospgrillage_obj, result_obj, component="Fy",
                    member=m, loadcase=loadcase, **kwargs,
                )
            )
        except (ValueError, KeyError, IndexError):
            pass
    return figs


def plot_def(ospgrillage_obj, result_obj=None, member=None, loadcase=None,
             backend="matplotlib", **kwargs):
    """
    Plot vertical deflection (y-displacement) for one or all main beams.

    When *member* is ``None``, iterates over the three main-beam member groups
    and returns a list of figures.

    :param ospgrillage_obj: Grillage model object.
    :param result_obj: xarray DataSet of results.
    :param member: Member name. If ``None``, plots all main beams.
    :param loadcase: Load case name. If ``None``, uses the first load case.
    :param backend: ``"matplotlib"`` (default, static) or ``"plotly"`` (interactive 3D).
    :param \\**kwargs: Forwarded to the underlying renderer.  See
        :func:`plot_defo` for accepted keyword arguments.
    :returns: Single axes when *member* is given, else list of axes.
        For ``backend="plotly"``, returns a single :class:`plotly.graph_objects.Figure`.
    """
    members = [member] if member else _MAIN_BEAM_MEMBERS
    if backend == "plotly":
        plotly_kw = {k: v for k, v in kwargs.items()
                     if k in ("figsize", "scale", "title")}
        plotly_kw["fig"] = kwargs.get("ax", None)
        return _plotly_3d_defo(
            ospgrillage_obj, result_obj, members=members,
            component="y", loadcase=loadcase, **plotly_kw,
        )
    # matplotlib path — filter to defo-compatible kwargs (no fill/alpha)
    defo_kw = {k: v for k, v in kwargs.items()
               if k in ("figsize", "ax", "scale", "title", "color", "show")}
    if member is not None:
        return plot_defo(
            ospgrillage_obj, result_obj, member=member,
            component="y", loadcase=loadcase, **defo_kw,
        )
    figs = []
    for m in _MAIN_BEAM_MEMBERS:
        try:
            figs.append(
                plot_defo(
                    ospgrillage_obj, result_obj, member=m,
                    component="y", loadcase=loadcase, **defo_kw,
                )
            )
        except (ValueError, KeyError, IndexError):
            pass
    return figs


# ---------------------------------------------------------------------------
# Model geometry plotting
# ---------------------------------------------------------------------------

# Colour palette for member groups (one per member type).
_MEMBER_COLOURS = {
    "edge_beam": "grey",
    "exterior_main_beam_1": "blue",
    "interior_main_beam": "green",
    "exterior_main_beam_2": "blue",
    "start_edge": "red",
    "end_edge": "red",
    "transverse_slab": "orange",
}
_DEFAULT_COLOUR = "black"


def _extract_mesh_data(grillage_obj):
    """Return element geometry grouped by member name.

    Returns
    -------
    dict[str, list[tuple]]
        ``{member_name: [(node_i_xyz, node_j_xyz, ele_tag), ...]}``
    all_node_coords : dict[int, list]
        ``{node_tag: [x, y, z]}``
    """
    mesh = grillage_obj.Mesh_obj
    node_spec = mesh.node_spec
    z_group_map = grillage_obj.common_grillage_element_z_group

    data = {}  # member_name -> list of (coord_i, coord_j, ele_tag)
    all_nodes = {}  # node_tag -> [x, y, z]

    def _add_elements(member_name, ele_list):
        entries = data.setdefault(member_name, [])
        for ele in ele_list:
            tag, ni, nj = ele[0], ele[1], ele[2]
            ci = node_spec[ni]["coordinate"]
            cj = node_spec[nj]["coordinate"]
            entries.append((ci, cj, tag))
            all_nodes[ni] = ci
            all_nodes[nj] = cj

    # Longitudinal members (edge_beam, exterior_main_beam_1, interior, exterior_2)
    for member_name in ["edge_beam", "exterior_main_beam_1",
                        "interior_main_beam", "exterior_main_beam_2"]:
        if member_name not in z_group_map:
            continue
        for z_idx in z_group_map[member_name]:
            if z_idx in mesh.z_group_to_ele:
                _add_elements(member_name, mesh.z_group_to_ele[z_idx])

    # Start / end edges
    for member_name in ["start_edge", "end_edge"]:
        if member_name not in z_group_map:
            continue
        for edge_idx in z_group_map[member_name]:
            if edge_idx in mesh.edge_group_to_ele:
                _add_elements(member_name, mesh.edge_group_to_ele[edge_idx])

    # Transverse slab
    _add_elements("transverse_slab", mesh.trans_ele)

    return data, all_nodes


def _plot_model_matplotlib(grillage_obj, *, figsize=None, ax=None,
                           title=_AUTO, show_nodes=True,
                           show_node_labels=False, show_element_labels=False,
                           color_by_member=True, show=False):
    """Render a 2-D plan view (x vs z) of the grillage mesh."""
    data, all_nodes = _extract_mesh_data(grillage_obj)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Draw elements grouped by member
    for member_name, elements in data.items():
        colour = _MEMBER_COLOURS.get(member_name, _DEFAULT_COLOUR) if color_by_member else "k"
        for ci, cj, etag in elements:
            ax.plot([ci[0], cj[0]], [ci[2], cj[2]], "-", color=colour, linewidth=1)
            if show_element_labels:
                mx = 0.5 * (ci[0] + cj[0])
                mz = 0.5 * (ci[2] + cj[2])
                ax.text(mx, mz, str(etag), fontsize=6, color="red",
                        ha="center", va="center")

    # Draw nodes
    if show_nodes or show_node_labels:
        for ntag, coord in all_nodes.items():
            if show_nodes:
                ax.plot(coord[0], coord[2], ".", color="k", markersize=3)
            if show_node_labels:
                ax.text(coord[0], coord[2], f" {ntag}", fontsize=5,
                        color="blue", ha="left", va="bottom")

    # Add legend (one entry per member group)
    from matplotlib.lines import Line2D
    handles = []
    for member_name in data:
        colour = _MEMBER_COLOURS.get(member_name, _DEFAULT_COLOUR) if color_by_member else "k"
        handles.append(Line2D([0], [0], color=colour, linewidth=1, label=member_name))
    if handles and color_by_member:
        ax.legend(handles=handles, fontsize=7, loc="best")

    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_aspect("equal")

    if title is _AUTO:
        ax.set_title("Grillage Model")
    elif title is not None:
        ax.set_title(title)

    fig.tight_layout()
    if show:
        plt.show()

    return ax


def _plot_model_plotly(grillage_obj, *, fig=None, figsize=None,
                       title=_AUTO, show_nodes=True,
                       show_node_labels=False, show_element_labels=False,
                       color_by_member=True, show=False):
    """Render an interactive 3-D model view with Plotly."""
    go = _import_plotly()
    if fig is None:
        fig = go.Figure()

    data, all_nodes = _extract_mesh_data(grillage_obj)

    colours = ["grey", "blue", "green", "blue", "red", "red", "orange"]
    colour_map = dict(zip(_MEMBER_COLOURS.keys(), colours))

    # Draw elements grouped by member
    for member_name, elements in data.items():
        colour = colour_map.get(member_name, "black") if color_by_member else "black"
        xs, ys, zs = [], [], []
        for ci, cj, etag in elements:
            xs.extend([ci[0], cj[0], None])
            ys.extend([ci[2], cj[2], None])
            zs.extend([ci[1], cj[1], None])
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs,
            mode="lines",
            line=dict(color=colour, width=3),
            name=member_name,
        ))

    # Node markers
    if show_nodes:
        ntags = list(all_nodes.keys())
        nx = [all_nodes[n][0] for n in ntags]
        ny = [all_nodes[n][2] for n in ntags]
        nz = [all_nodes[n][1] for n in ntags]
        text = [str(n) for n in ntags] if show_node_labels else None
        mode = "markers+text" if show_node_labels else "markers"
        fig.add_trace(go.Scatter3d(
            x=nx, y=ny, z=nz,
            mode=mode,
            marker=dict(size=2, color="black"),
            text=text,
            textposition="top center",
            textfont=dict(size=8),
            name="nodes",
            showlegend=False,
        ))

    # Element labels at midpoints
    if show_element_labels:
        ex, ey, ez, etexts = [], [], [], []
        for member_name, elements in data.items():
            for ci, cj, etag in elements:
                ex.append(0.5 * (ci[0] + cj[0]))
                ey.append(0.5 * (ci[2] + cj[2]))
                ez.append(0.5 * (ci[1] + cj[1]))
                etexts.append(str(etag))
        fig.add_trace(go.Scatter3d(
            x=ex, y=ey, z=ez,
            mode="text",
            text=etexts,
            textfont=dict(size=7, color="red"),
            showlegend=False,
        ))

    # Layout
    layout_kw = dict(
        scene=dict(
            xaxis_title="x (m)",
            yaxis_title="z (m)",
            zaxis_title="y (m)",
            aspectmode="data",
        ),
        legend_title="Member",
    )
    if title is _AUTO:
        layout_kw["title"] = "Grillage Model"
    elif title is not None:
        layout_kw["title"] = title
    if figsize is not None:
        layout_kw["width"] = figsize[0] * 100
        layout_kw["height"] = figsize[1] * 100

    fig.update_layout(**layout_kw)

    if show:
        fig.show()

    return fig


def plot_model(grillage_obj, *, backend="matplotlib", **kwargs):
    """Plot the grillage mesh geometry.

    :param grillage_obj: Grillage model (must have been created with
        :meth:`~ospgrillage.osp_grillage.OspGrillage.create_osp_model`).
    :type grillage_obj: OspGrillage
    :param backend: ``"matplotlib"`` (default, 2-D plan view) or
        ``"plotly"`` (interactive 3-D).
    :type backend: str
    :param figsize: Figure size in inches ``(width, height)``.
        Ignored when *ax* is provided.
    :type figsize: tuple, optional
    :param ax: Existing matplotlib Axes to plot on (matplotlib backend only).
    :type ax: :class:`~matplotlib.axes.Axes`, optional
    :param fig: Existing Plotly Figure to add traces to (plotly backend only).
    :type fig: :class:`plotly.graph_objects.Figure`, optional
    :param title: Plot title.  Default auto-generates ``"Grillage Model"``.
        Pass a string to override, or ``None`` to suppress.
    :type title: str or None
    :param show_nodes: Show node markers.  Default ``True``.
    :type show_nodes: bool
    :param show_node_labels: Annotate node tags.  Default ``False``.
    :type show_node_labels: bool
    :param show_element_labels: Annotate element tags.  Default ``False``.
    :type show_element_labels: bool
    :param color_by_member: Colour elements by member group.  Default ``True``.
    :type color_by_member: bool
    :param show: If ``True``, call ``plt.show()`` before returning.
    :type show: bool
    :returns: Matplotlib axes (use ``ax.get_figure()`` for the figure) or
        Plotly figure.
    :rtype: :class:`~matplotlib.axes.Axes` or
        :class:`plotly.graph_objects.Figure`
    """
    if backend == "plotly":
        plotly_kw = {k: v for k, v in kwargs.items()
                     if k in ("fig", "figsize", "title", "show_nodes",
                              "show_node_labels", "show_element_labels",
                              "color_by_member", "show")}
        # Map ax → fig for consistency with other convenience wrappers
        if "ax" in kwargs and "fig" not in plotly_kw:
            plotly_kw["fig"] = kwargs["ax"]
        return _plot_model_plotly(grillage_obj, **plotly_kw)
    elif backend == "matplotlib":
        mpl_kw = {k: v for k, v in kwargs.items()
                  if k in ("figsize", "ax", "title", "show_nodes",
                           "show_node_labels", "show_element_labels",
                           "color_by_member", "show")}
        return _plot_model_matplotlib(grillage_obj, **mpl_kw)
    else:
        raise ValueError(
            f"Unknown backend: {backend!r}. Use 'matplotlib' or 'plotly'."
        )


class PostProcessor:
    """Class to post-process the results from an of ospgrillage analysis result.

    As of version 0.2.0, ospgrillage compiles results using Xarray module.
    """

    def __init__(self, grillage, result):  # Union[xr.DataArray, xr.Dataset]
        """Init the class"""
        # store main vars
        self.grillage = grillage
        self.result = result

        # init vars
        self.shape_function_obj = ShapeFunction()

    def get_arbitrary_displacements(
        self, point: list, shape_function_type: str = "linear", component: str = "y"
    ):
        """Returns displacement (translation and rotational) from an arbitrary point by interpolation

        param point: list of coordinate. Default three elements [x,y=0,z]
        type point: list
        param component: Displacement component. Default "y"
        type component: str
        param shape_function_type: The shape function for interpolation. Default "linear"
        type shape_function_type: str
        """
        node_displacements = []
        node_coordinate = []

        # get the list of four nodes where arbitrary point lies
        nodes, grid_number = self.grillage._get_point_load_nodes(
            point=point
        )  # list of nodes
        if nodes is None:
            raise ValueError("Point is outside bridge mesh")

        # get results of each node of four nodes
        for node in nodes:
            node_displacements.append(
                self.result.sel(
                    Component=component,
                    Node=node,
                )
                .displacements.to_numpy()
                .tolist()[0]
            )
            node_coordinate.append(self.grillage.get_nodes(number=node))

        # interpolate for vertical displacement
        x = np.array([coord[0] for coord in node_coordinate])
        z = np.array([coord[2] for coord in node_coordinate])

        # get natural coordinate of point in grid
        eta, zeta = solve_zeta_eta(
            xp=point[0],
            zp=point[2],
            x1=x[0],
            z1=z[0],
            x2=x[1],
            z2=z[1],
            x3=x[2],
            z3=z[2],
            x4=x[3],
            z4=z[3],
        )
        if shape_function_type == "linear":
            shape_func = self.shape_function_obj.linear_shape_function(
                eta=eta, zeta=zeta
            )
        else:
            shape_func, _, _ = self.shape_function_obj.hermite_shape_function_2d(
                eta=eta, zeta=zeta
            )

        return sum([a * b for a, b, in zip(shape_func, node_displacements)])
