# SUPPORT SUMMARY
# Maintain this section whenever functionality, supported elements/analyses,
# stress formats, extraction modes, or limitations change.
#
# Functionality:
# - ResponseExtractor.extract_von_mises(model_part, mode="average") extracts
#   one von Mises stress per element as a NumPy array in element iteration order.
# - ResponseExtractor.extract_signed_von_mises(model_part, mode="average")
#   returns the von Mises magnitude multiplied by the sign of the principal
#   stress with the largest absolute value. Exact opposite-sign magnitude ties
#   use the positive principal stress. A zero von Mises magnitude remains zero.
# - ResponseExtractor.extract_tresca(model_part, mode="average") extracts one
#   Tresca equivalent stress per element: max(principal stresses) minus min.
#   This is twice the maximum shear stress, with the same output array format.
# - ResponseExtractor.extract_maximum_shear_stress(model_part, mode="average")
#   returns (sigma1 - sigma3)/2, where sigma1 >= sigma2 >= sigma3.
# - ResponseExtractor.extract_signed_maximum_shear_stress(model_part,
#   mode="average") gives maximum shear the sign of the principal stress with
#   the largest absolute value, using the same tie rule as signed von Mises.
# - ResponseExtractor.extract_principle_stresses(model_part, mode="average")
#   returns a list of [sigma1, sigma2, sigma3] lists in element iteration order,
#   with sigma1 >= sigma2 >= sigma3 (descending algebraic principal values).
# - ResponseExtractor.extract_principle_strains(model_part, mode="average")
#   returns [epsilon1, epsilon2, epsilon3] lists, with epsilon1 >= epsilon2 >=
#   epsilon3. Strain outputs and additional restrictions are described below.
# - ResponseExtractor.extract_von_mises_strain(model_part, mode="average")
#   applies the Abaqus/CAE MISES tensor invariant to small strain:
#   sqrt(3/2 * epsilon_dev:epsilon_dev). This is not accumulated equivalent
#   plastic strain (PEEQ) or the sqrt(2/3 * epsilon_dev:epsilon_dev) definition.
# - Uses the Abaqus invariant sqrt(3/2 s:s), where s is the deviatoric Cauchy
#   stress tensor, for von Mises. Stress and strain formats and their support
#   limits are documented below. Call after computing the element's
#   displacement/stress state; this helper does not solve the analysis.
#
# Supported analyses and elements (including first AND second order):
# - 2D plane stress: continuum triangles and quadrilaterals with a constitutive
#   law declaring PLANE_STRESS_LAW. Examples: SmallDisplacementElement2D3N,
#   2D4N (first order), and 2D6N, 2D8N, 2D9N (second order).
# - 2D plane strain: the same continuum shapes and interpolation orders, with
#   a constitutive law declaring PLANE_STRAIN_LAW. Examples include the same
#   SmallDisplacementElement2D3N/4N/6N/8N/9N elements, subject to sigma_zz below.
# - 3D solids: tetrahedra, pyramids, prisms, and hexahedra, with a constitutive
#   law having working dimension 3 and strain size 6. Examples include
#   SmallDisplacementElement3D4N/5N/6N/8N (first order) and
#   SmallDisplacementElement3D10N/13N/15N/20N/27N (second order).
# - Support is determined by geometry, law flags, and available Cauchy stress
#   output, not a fixed element-name list. Complete stress output is not
#   restricted to linear-elastic analyses; reduced plane strain is restricted.
#
# Stress output and plane-strain limitations:
# - Requests CAUCHY_STRESS_TENSOR for plane stress/3D solids and
#   CAUCHY_STRESS_VECTOR for plane strain, preserving available sigma_zz.
# - Accepts full 3x3 tensors and six-component Kratos Voigt vectors ordered
#   [xx, yy, zz, xy, yz, xz]. In 2D also accepts 2x2 tensors, three-component
#   [xx, yy, xy], and four-component [xx, yy, zz, xy] vectors.
# - Full supplied tensors retain every component. Reduced plane-stress output
#   has zero out-of-plane stresses by the plane-stress assumption.
# - Plane strain requires sigma_zz. If omitted, sigma_zz = nu*(sigma_xx+sigma_yy)
#   is recovered ONLY for SmallDisplacementElement... or
#   SmallDisplacementBbarElement... with exactly LinearElasticPlaneStrain2DLaw
#   and no initial stress or strain. Other missing-sigma_zz cases are unsupported.
#
# Strain output and restrictions:
# - Strain extractors currently support SmallDisplacementElement... and 3D
#   SmallDisplacementBbarElement... with the built-in LinearElastic3DLaw,
#   LinearElasticPlaneStress2DLaw, or LinearElasticPlaneStrain2DLaw only.
# - Read GREEN_LAGRANGE_STRAIN_VECTOR, which these formulations return using
#   infinitesimal kinematics. Finite-strain formulations are unsupported.
# - Voigt strain order is [xx, yy, zz, gamma_xy, gamma_yz, gamma_xz], with
#   engineering shear divided by TWO when constructing the strain tensor.
#   The 2D orders are [xx, yy, gamma_xy] and [xx, yy, zz, gamma_xy]. Full
#   tensors already contain physical shear strain and are used unchanged.
# - Reduced plane strain has epsilon_zz = 0 by its kinematic constraint.
# - Reduced plane stress uses epsilon_zz = 0 to match Abaqus/CAE strain FIELD
#   output. Physical thickness strain may be nonzero; this is the stored-field
#   convention, not a reconstruction of thickness contraction. Explicit
#   out-of-plane strain is preserved when present in full supplied output.
#
# Extraction modes (available for all extractors):
# - average (default): computes the selected response at each GP, then averages
#   the scalars or each ordered principal value; does not average tensors first.
#   Signed responses average signed GP values, allowing cancellation.
# - Gp_max for unsigned scalar responses: returns the largest GP value.
# - Gp_max for signed responses: use the GP with the largest response magnitude
#   and retain its sign. Ties use the first GP. For example, -100 wins over +20.
# - Gp_max for extract_principle_stresses: use the Gauss point with the largest
#   sigma1 in each element. Return [sigma1, sigma2, sigma3] from that SAME Gauss
#   point; sigma2 and sigma3 are not maximized independently. Ties use the first GP.
# - Gp_max for extract_principle_strains: use the GP with the largest epsilon1
#   and return all three principal strains from that SAME GP. Ties use the first GP.
# - center_gp: reads the zero-based index from center_gauss_point_lookup.py,
#   matching a listed element/geometry and its ACTIVE GP count to GI_GAUSS_1..5.
#   Uses the mapped reference-element centroid defined in the source catalogue;
#   no coordinate comparisons, tolerance, or nearest-point calculation is used.
#   Unknown elements/rules raise NotImplementedError. The old coordinate-based
#   implementation is retained below as commented reference code.
#   Without a center GP, warns and calls the center_element fallback.
# - center_element: NOT IMPLEMENTED; raises NotImplementedError, including
#   when reached through the center_gp fallback.
#
# Unsupported: trusses, beams, shells, membranes, axisymmetric analyses, and
# other element/law combinations failing the support checks. These raise
# NotImplementedError; incompatible stress-array shapes raise ValueError.

from pathlib import Path
from runpy import run_path

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np


# Load the sibling data file once, including when responses.py is loaded directly.
_CENTER_GAUSS_POINT_LOOKUP = run_path(
    str(Path(__file__).with_name("center_gauss_point_lookup.py"))
)["CENTER_GAUSS_POINT_LOOKUP"]


class ResponseExtractor:

    def extract_von_mises(self, model_part, mode="average"):
        """Return a NumPy array of von Mises stresses in element iteration order.

        ``average`` averages the von Mises values at all Gauss points.
        ``Gp_max`` returns the maximum von Mises value among all Gauss points.
        ``center_gp`` uses a Gauss point at the geometry center, falling back
        to ``center_element`` with a warning if none exists.
        ``center_element`` is reserved and not implemented yet.
        Supports plane-stress/plane-strain continua and 3D solids. Uses the
        Abaqus Mises invariant sqrt(3/2 s:s), including sigma_zz in plane strain.
        Reduced plane-strain output can only be recovered for the built-in
        small-strain isotropic elastic law without initial stress or strain.
        """
        return self._extract_response(model_part, mode, self._von_mises)

    def extract_signed_von_mises(self, model_part, mode="average"):
        """Return signed von Mises stresses in element iteration order.

        The sign comes from the principal Cauchy stress with the largest
        absolute value; exact opposite-sign ties use the positive principal.
        ``average`` averages the signed Gauss-point values.
        ``Gp_max`` returns the value with the largest von Mises magnitude,
        retaining its sign; ties use the first GP.
        ``center_gp`` uses a center GP or warns and falls back to the unfinished
        ``center_element`` mode, which raises NotImplementedError.
        Returns a NumPy array with the same supported elements, analyses,
        stress formats, and plane-strain sigma_zz recovery as extract_von_mises.
        """
        return self._extract_response(
            model_part, mode, self._signed_von_mises, gp_max_absolute=True
        )

    def extract_tresca(self, model_part, mode="average"):
        """Return a NumPy array of Tresca stresses in element iteration order.

        Tresca equivalent stress is the largest minus the smallest principal
        Cauchy stress (twice the maximum shear stress).
        ``average`` averages the Tresca values at all Gauss points.
        ``Gp_max`` returns the maximum Tresca value among all Gauss points.
        ``center_gp`` uses a center Gauss point or warns and falls back to
        ``center_element``, which is not implemented yet.
        Supports the same elements, analyses, stress formats, and plane-strain
        sigma_zz recovery as extract_von_mises.
        """
        return self._extract_response(model_part, mode, self._tresca)

    def extract_maximum_shear_stress(self, model_part, mode="average"):
        """Return maximum shear stress (sigma1 - sigma3)/2 per element.

        ``average`` averages GP maximum-shear values; ``Gp_max`` takes the
        largest GP value. ``center_gp`` uses a center GP or warns and falls
        back to the unimplemented ``center_element`` mode.
        Uses the same supported elements, analyses, and stress formats as
        extract_von_mises, returning a NumPy array in element iteration order.
        """
        return self._extract_response(model_part, mode, self._maximum_shear_stress)

    def extract_signed_maximum_shear_stress(self, model_part, mode="average"):
        """Return signed maximum shear stress in element iteration order.

        Uses the sign of the largest-magnitude principal Cauchy stress; exact
        opposite-sign ties favor tension. ``average`` averages signed GP
        values; ``Gp_max`` selects the largest magnitude and retains its sign.
        ``center_gp`` and the unfinished ``center_element`` mode behave as in
        extract_von_mises. Returns a NumPy array with the same support limits.
        """
        return self._extract_response(
            model_part, mode, self._signed_maximum_shear_stress, gp_max_absolute=True
        )

    def extract_principle_stresses(self, model_part, mode="average"):
        """Return one [sigma1, sigma2, sigma3] list per element, largest first.

        ``average`` averages each ordered principal value across Gauss points.
        ``Gp_max`` returns the complete triplet from the GP with the largest
        sigma1; ties use the first GP. Values are ordered algebraically as
        sigma1 >= sigma2 >= sigma3, including compressive and out-of-plane stress.
        ``center_gp`` uses a center GP or warns and falls back to the unfinished
        ``center_element`` mode, which raises NotImplementedError.
        Uses the same supported elements, analyses, and stress formats as
        extract_von_mises. The outer list follows element iteration order.
        """
        return self._extract_response(
            model_part, mode, self._principal_stresses
        ).tolist()

    def extract_principle_strains(self, model_part, mode="average"):
        """Return [epsilon1, epsilon2, epsilon3] lists in element iteration order.

        Values are sorted epsilon1 >= epsilon2 >= epsilon3. ``average``
        averages each ordered GP principal value. ``Gp_max`` returns the full
        triplet at the GP with the largest epsilon1; ties use the first GP.
        ``center_gp`` uses a center GP or warns and calls the unfinished
        ``center_element`` fallback. Strain support is restricted to the
        small-displacement linear-elastic formulations in the opening summary.
        Missing epsilon_zz is zero for reduced plane strain by kinematics and
        for reduced plane stress by the Abaqus/CAE field-output convention.
        Explicitly supplied out-of-plane strain is preserved.
        """
        return self._extract_response(
            model_part, mode, self._principal_strains, strain_response=True
        ).tolist()

    def extract_von_mises_strain(self, model_part, mode="average"):
        """Return the Abaqus/CAE MISES tensor invariant of small strain.

        Computes sqrt(3/2 * epsilon_dev:epsilon_dev) after converting engineering
        shear to physical tensor shear. This is a tensor invariant, not PEEQ.
        ``average`` averages the GP invariants; ``Gp_max`` takes the largest.
        ``center_gp`` uses a center GP or warns and calls the unimplemented
        ``center_element`` fallback. Returns a NumPy array in element order.
        Uses the strain support and out-of-plane field convention at the top.

        Abaqus FieldOutput documents MISES and engineering-tensor handling:
        https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKERRefMap/simaker-c-fieldoutputpyc.htm
        """
        return self._extract_response(
            model_part, mode, self._von_mises_strain, strain_response=True
        )

    def _extract_response(self, model_part, mode, response_measure,
                          gp_max_absolute=False, strain_response=False):
        """Extract scalar or principal stress/strain values using the GP mode."""
        if mode not in ("center_gp", "average", "center_element", "Gp_max"):
            raise ValueError(
                "mode must be 'center_gp', 'average', 'center_element', or 'Gp_max'."
            )
        values = []
        for element in model_part.Elements:
            analysis_type = self._get_analysis_type(element)
            if strain_response:
                self._check_strain_support(element)
            if mode == "center_element":
                values.append(self._extract_center_element(element))
                continue
            if mode == "center_gp":
                gp_index = self._get_center_gp_index(element)
                if gp_index is None:
                    KM.Logger.PrintWarning(
                        "ResponseExtractor",
                        f"Element {element.Id} has no Gauss point at its center; "
                        "falling back to 'center_element'."
                    )
                    values.append(self._extract_center_element(element))
                    continue
            if strain_response:
                variable = KM.GREEN_LAGRANGE_STRAIN_VECTOR
            else:
                # A plane-strain vector can include sigma_zz even when the
                # element's 2D tensor output cannot store it.
                variable = (KM.CAUCHY_STRESS_VECTOR if analysis_type == "plane_strain"
                            else KM.CAUCHY_STRESS_TENSOR)
            gp_data = element.CalculateOnIntegrationPoints(variable, model_part.ProcessInfo)
            if not gp_data:
                raise RuntimeError(
                    f"Element {element.Id} returned no Gauss-point "
                    f"{'strain' if strain_response else 'stress'} tensors."
                )
            if mode == "center_gp" and len(gp_data) != len(element.GetIntegrationPoints()):
                raise RuntimeError(
                    f"Element {element.Id}: Gauss-point response count does not "
                    "match the active integration-point count."
                )
            poisson_ratio = None
            reduced_output = any(np.asarray(data).shape in ((3,), (2, 2)) for data in gp_data)
            if reduced_output and not strain_response and analysis_type == "plane_strain":
                poisson_ratio = self._plane_strain_poisson_ratio(element, model_part.ProcessInfo)
            measure_args = (analysis_type,) if strain_response else (analysis_type, poisson_ratio)
            if mode in ("average", "Gp_max"):
                gp_values = np.asarray([
                    response_measure(data, *measure_args)
                    for data in gp_data
                ])
                if mode == "average":
                    value = np.mean(gp_values, axis=0)
                elif gp_values.ndim == 2:
                    value = gp_values[np.argmax(gp_values[:, 0])]
                elif gp_max_absolute:
                    value = gp_values[np.argmax(np.abs(gp_values))]
                else:
                    value = np.max(gp_values)
            else:
                value = response_measure(gp_data[gp_index], *measure_args)
            values.append(value)
        return np.asarray(values)

    @staticmethod
    def _check_strain_support(element):
        """Limit strain measures to the supported small-strain elastic cases."""
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(element)
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        supported_formulation = (name.startswith("SmallDisplacementElement")
                                 or name.startswith("SmallDisplacementBbarElement3D"))
        if (supported_formulation
                and type(law) in (SMA.LinearElastic3DLaw, SMA.LinearElasticPlaneStress2DLaw,
                                  SMA.LinearElasticPlaneStrain2DLaw)):
            return
        raise NotImplementedError(
            f"Element {element.Id} ({name}): strain extraction currently requires "
            "a small-displacement continuum with a built-in linear-elastic law."
        )

    @staticmethod
    def _get_analysis_type(element):
        """Identify supported continuum analyses from geometry and law flags."""
        geometry = element.GetGeometry()
        dimension = geometry.LocalSpaceDimension()
        if dimension in (2, 3) and dimension == geometry.WorkingSpaceDimension():
            if element.Properties.Has(KM.CONSTITUTIVE_LAW):
                law = element.Properties[KM.CONSTITUTIVE_LAW]
                features = KM.ConstitutiveLawFeatures()
                law.GetLawFeatures(features)
                options = features.GetOptions()
                if dimension == 2:
                    if options.Is(KM.ConstitutiveLaw.PLANE_STRESS_LAW):
                        return "plane_stress"
                    if options.Is(KM.ConstitutiveLaw.PLANE_STRAIN_LAW):
                        return "plane_strain"
                elif law.WorkingSpaceDimension() == 3 and law.GetStrainSize() == 6:
                    return "solid"
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(element)
        raise NotImplementedError(
            f"Element {element.Id} ({name}): unsupported element or analysis; "
            "only plane stress, plane strain and 3D solids are supported."
        )

    @staticmethod
    def _plane_strain_poisson_ratio(element, process_info):
        """Allow sigma_zz = nu (sigma_xx + sigma_yy) only where it is valid."""
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(element)
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        if (name.startswith(("SmallDisplacementElement", "SmallDisplacementBbarElement"))
                and type(law) is SMA.LinearElasticPlaneStrain2DLaw):
            for variable in (KM.INITIAL_STRESS_VECTOR, KM.INITIAL_STRAIN_VECTOR):
                values = element.CalculateOnIntegrationPoints(variable, process_info)
                if not values or any(np.any(np.asarray(value) != 0.0) for value in values):
                    break
            else:
                return element.Properties[KM.POISSON_RATIO]
        raise NotImplementedError(
            f"Element {element.Id} ({name}): plane-strain output is missing sigma_zz. "
            "Recovery requires the built-in small-strain isotropic elastic law "
            "without initial stress or strain."
        )

    @staticmethod
    def _get_center_gp_index(element):
        """Read the zero-based center index from the standard-rule catalogue."""
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(element)
        geometry_type = element.GetGeometry().GetGeometryType().name.removeprefix("Kratos_")
        element_data = _CENTER_GAUSS_POINT_LOOKUP["elements"].get(geometry_type)
        if element_data is None or name not in element_data["element_names"]:
            raise NotImplementedError(
                f"Element {element.Id} ({name}, {geometry_type}): no center-GP "
                "lookup for this element. Only the listed standard geometry-based "
                "GI_GAUSS_1 through GI_GAUSS_5 formulations are supported."
            )
        number_of_points = len(element.GetIntegrationPoints())
        rules = _CENTER_GAUSS_POINT_LOOKUP["families"][element_data["family"]]["rules_by_gp_count"]
        rule = rules.get(number_of_points)
        if rule is None:
            raise NotImplementedError(
                f"Element {element.Id} ({name}): {number_of_points} integration "
                "points do not match a catalogued GI_GAUSS rule for this geometry."
            )
        return rule["center_gp_index"]

    # OLD IMPLEMENTATION: coordinate comparison with a numerical tolerance.
    # Preserved for reference; center_gp now uses the catalogue above.
    # @staticmethod
    # def _get_center_gp_index(element):
    #     """Find the center GP using the element's actual integration rule."""
    #     # GetIntegrationPoints uses the element geometry and active integration
    #     # method, accounting for the element type and selected integration order.
    #     points = element.GetIntegrationPoints()
    #     if not points:
    #         return None
    #     geometry = element.GetGeometry()
    #     center = geometry.Center()
    #     center_coordinates = np.array([center.X, center.Y, center.Z])
    #     tolerance = 1e-10 * geometry.Length()
    #     middle = len(points) // 2
    #     # Odd tensor-product rules usually store the center at the middle index.
    #     # Other rules can store it elsewhere, so verify coordinates in both cases.
    #     indices = [middle] + [i for i in range(len(points)) if i != middle]
    #     for gp_index in indices:
    #         distance = np.linalg.norm(np.asarray(points[gp_index]) - center_coordinates)
    #         if distance <= tolerance:
    #             return gp_index
    #     return None

    def _extract_center_element(self, element):
        # Reserved for element-center stress recovery.
        raise NotImplementedError("'center_element' is not implemented yet.")

    @staticmethod
    def _von_mises(stress, analysis_type="solid", poisson_ratio=None):
        """Compute the Abaqus Mises invariant from tensor or Kratos Voigt output.

        See Abaqus Conventions, Stress Invariants:
        https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMODRefMap/simamod-c-conventions.htm
        """
        tensor = ResponseExtractor._stress_tensor(stress, analysis_type, poisson_ratio)
        deviatoric_stress = tensor - np.trace(tensor) / 3.0 * np.eye(3)
        return np.sqrt(1.5 * np.sum(deviatoric_stress**2))

    @staticmethod
    def _signed_von_mises(stress, analysis_type="solid", poisson_ratio=None):
        """Use the sign of the largest-magnitude principal stress; ties favor tension."""
        principal_stresses = ResponseExtractor._principal_stresses(
            stress, analysis_type, poisson_ratio
        )
        dominant_stress = principal_stresses[np.argmax(np.abs(principal_stresses))]
        return np.sign(dominant_stress) * ResponseExtractor._von_mises(
            stress, analysis_type, poisson_ratio
        )

    @staticmethod
    def _tresca(stress, analysis_type="solid", poisson_ratio=None):
        """Return the largest difference between principal Cauchy stresses."""
        tensor = ResponseExtractor._stress_tensor(stress, analysis_type, poisson_ratio)
        principal_stresses = np.linalg.eigvalsh(tensor)
        return principal_stresses[-1] - principal_stresses[0]

    @staticmethod
    def _maximum_shear_stress(stress, analysis_type="solid", poisson_ratio=None):
        """Return half the largest principal-stress difference."""
        return 0.5 * ResponseExtractor._tresca(stress, analysis_type, poisson_ratio)

    @staticmethod
    def _signed_maximum_shear_stress(stress, analysis_type="solid", poisson_ratio=None):
        """Give maximum shear the sign of the largest-magnitude principal stress."""
        principal_stresses = ResponseExtractor._principal_stresses(
            stress, analysis_type, poisson_ratio
        )
        dominant_stress = principal_stresses[np.argmax(np.abs(principal_stresses))]
        return np.sign(dominant_stress) * 0.5 * (
            principal_stresses[0] - principal_stresses[-1]
        )

    @staticmethod
    def _principal_stresses(stress, analysis_type="solid", poisson_ratio=None):
        """Return principal Cauchy stresses ordered sigma1 >= sigma2 >= sigma3."""
        tensor = ResponseExtractor._stress_tensor(stress, analysis_type, poisson_ratio)
        return np.linalg.eigvalsh(tensor)[::-1]

    @staticmethod
    def _principal_strains(strain, analysis_type="solid"):
        """Return ordered small-strain principal values, including epsilon_zz."""
        tensor = ResponseExtractor._strain_tensor(strain, analysis_type)
        return np.linalg.eigvalsh(tensor)[::-1]

    @staticmethod
    def _von_mises_strain(strain, analysis_type="solid"):
        """Apply the Abaqus MISES field invariant to the physical strain tensor."""
        tensor = ResponseExtractor._strain_tensor(strain, analysis_type)
        deviatoric_strain = tensor - np.trace(tensor) / 3.0 * np.eye(3)
        return np.sqrt(1.5 * np.sum(deviatoric_strain**2))

    @staticmethod
    def _strain_tensor(strain, analysis_type="solid"):
        """Convert strain to 3D using Abaqus/CAE reduced plane-stress fields.

        https://docs.software.vt.edu/abaqusv2025/English/SIMACAEOUTRefMap/simaout-c-stdoutputvar.htm
        """
        strain = np.asarray(strain)
        if analysis_type not in ("plane_stress", "plane_strain", "solid"):
            raise NotImplementedError(f"Unsupported analysis type: {analysis_type}.")
        if strain.shape == (3, 3):
            return strain
        if strain.shape == (6,):
            xx, yy, zz, gamma_xy, gamma_yz, gamma_xz = strain
            return np.array([[xx, gamma_xy / 2.0, gamma_xz / 2.0],
                             [gamma_xy / 2.0, yy, gamma_yz / 2.0],
                             [gamma_xz / 2.0, gamma_yz / 2.0, zz]])
        if analysis_type in ("plane_stress", "plane_strain"):
            if strain.shape == (4,):
                xx, yy, zz, gamma_xy = strain
                return np.array([[xx, gamma_xy / 2.0, 0.0],
                                 [gamma_xy / 2.0, yy, 0.0], [0.0, 0.0, zz]])
            tensor = np.zeros((3, 3))
            if strain.shape == (2, 2):
                tensor[:2, :2] = strain
            elif strain.shape == (3,):
                xx, yy, gamma_xy = strain
                tensor[:2, :2] = [[xx, gamma_xy / 2.0], [gamma_xy / 2.0, yy]]
            else:
                raise ValueError(f"Unsupported {analysis_type} strain shape: {strain.shape}.")
            # Missing zz is zero for plane strain by kinematics and for plane
            # stress by the Abaqus/CAE strain-field convention (not physical
            # thickness-strain recovery).
            return tensor
        raise ValueError("A 3D solid requires a 3x3 strain tensor or six-component Voigt vector.")

    @staticmethod
    def _stress_tensor(stress, analysis_type, poisson_ratio):
        """Convert supported stress outputs to a complete 3D stress tensor."""
        stress = np.asarray(stress)
        if analysis_type not in ("plane_stress", "plane_strain", "solid"):
            raise NotImplementedError(f"Unsupported analysis type: {analysis_type}.")
        if stress.shape == (3, 3):
            tensor = stress
        elif stress.shape == (6,):
            # Kratos ordering differs from Abaqus/Standard: xx, yy, zz, xy, yz, xz.
            xx, yy, zz, xy, yz, xz = stress
            tensor = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        elif analysis_type in ("plane_stress", "plane_strain"):
            tensor = np.zeros((3, 3))
            if stress.shape == (4,):
                # Kratos ordering: xx, yy, zz, xy.
                xx, yy, zz, xy = stress
                tensor = np.array([[xx, xy, 0.0], [xy, yy, 0.0], [0.0, 0.0, zz]])
            else:
                if stress.shape == (2, 2):
                    tensor[:2, :2] = stress
                elif stress.shape == (3,):
                    xx, yy, xy = stress
                    tensor[:2, :2] = [[xx, xy], [xy, yy]]
                else:
                    raise ValueError(f"Unsupported {analysis_type} stress shape: {stress.shape}.")
                # In plane stress the omitted components are zero by definition.
                # Plane strain requires sigma_zz from the material response.
                if analysis_type == "plane_strain":
                    if poisson_ratio is None:
                        raise NotImplementedError("Plane-strain stress extraction requires sigma_zz.")
                    tensor[2, 2] = poisson_ratio * (tensor[0, 0] + tensor[1, 1])
        else:
            raise ValueError("A 3D solid requires a 3x3 tensor or six-component Voigt vector.")
        return tensor


