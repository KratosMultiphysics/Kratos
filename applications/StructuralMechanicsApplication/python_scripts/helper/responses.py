# SUPPORT SUMMARY
# Maintain this section whenever functionality, supported elements/analyses,
# stress formats, extraction modes, or limitations change.
#
# Functionality:
# - ResponseExtractor.extract_von_mises(model_part, mode="average") extracts
#   one von Mises stress per element as a NumPy array in element iteration
#   order.
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
#   returns [sigma1, sigma2, sigma3] lists in element iteration order,
#   with sigma1 >= sigma2 >= sigma3 (descending algebraic principal values).
# - ResponseExtractor.extract_principle_strains(model_part, mode="average")
#   returns [epsilon1, epsilon2, epsilon3] lists, with epsilon1 >= epsilon2 >=
#   epsilon3. Strain outputs and additional restrictions are described below.
# - ResponseExtractor.extract_von_mises_strain(model_part, mode="average")
#   applies the Abaqus/CAE MISES tensor invariant to small strain:
#   sqrt(3/2 * epsilon_dev:epsilon_dev). This is not accumulated equivalent
#   plastic strain (PEEQ) or sqrt(2/3 * epsilon_dev:epsilon_dev).
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
#   SmallDisplacementElement2D3N/4N/6N/8N/9N, subject to sigma_zz below.
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
# - Plane strain requires sigma_zz. If omitted, recover
#   sigma_zz = nu*(sigma_xx+sigma_yy) ONLY for SmallDisplacementElement... or
#   SmallDisplacementBbarElement... with exactly LinearElasticPlaneStrain2DLaw
#   and no initial stress or strain. Other missing-sigma_zz cases are
#   unsupported.
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
#   the scalars or each ordered principal value; does not average tensors
#   first.
#   Signed responses average signed GP values, allowing cancellation.
# - Gp_max for unsigned scalar responses: returns the largest GP value.
# - Gp_max for signed responses: use the GP with the largest response magnitude
#   and retain its sign. Ties use the first GP. For example, -100 wins over
#   +20.
# - Gp_max for extract_principle_stresses: use the Gauss point with the largest
#   sigma1 in each element. Return [sigma1, sigma2, sigma3] from that SAME
#   Gauss point; sigma2 and sigma3 are not maximized independently. Ties use
#   the first GP.
# - Gp_max for extract_principle_strains: use the GP with the largest epsilon1
#   and return all three principal strains from that SAME GP. Ties use the
#   first GP.
# - center_gp: reads the zero-based index from center_gauss_point_lookup.py,
#   matching a listed element/geometry and its ACTIVE GP count to
#   GI_GAUSS_1..5.
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
from typing import List, NoReturn, Optional

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np


# Load the sibling data file once, including for direct responses.py loading.
_CENTER_GAUSS_POINT_LOOKUP = run_path(
    str(Path(__file__).with_name("center_gauss_point_lookup.py"))
)["CENTER_GAUSS_POINT_LOOKUP"]


class ResponseExtractor:
    """Post-process continuum stresses and supported small-strain responses.

    Responses use the current element state; this class does not solve an
    analysis. Stress units follow the model's consistent unit system, and
    strains are dimensionless. Results follow ``model_part.Elements`` order.
    See the opening support summary for element, law, and output restrictions.
    """

    def extract_von_mises(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> np.ndarray:
        """Extract one von Mises Cauchy stress per continuum element.

        Mathematical Formulation
        ------------------------
        s = sigma - trace(sigma)/3 * I
        sigma_vm = sqrt(3/2 * sum_ij(s_ij**2))

        Here sigma is the Cauchy stress tensor, s its deviatoric part, and I
        the 3D identity. The invariant includes out-of-plane stress.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized elements and a computed solution state.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': arithmetic mean of GP von Mises values,
            without quadrature weights. 'Gp_max' takes the largest GP value.
            'center_gp' selects the mapped reference-element centroid; if its
            rule has no center GP, it warns and calls 'center_element'.
            'center_element' is reserved and raises when reached.

        Returns
        -------
        numpy.ndarray, shape (n_elements,)
            Nonnegative stresses in the model's stress units, in element
            iteration order. An empty model part returns an empty array.

        Raises
        ------
        ValueError
            If the mode or a stress-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If an element, law, analysis, or center rule is unsupported,
            sigma_zz cannot be recovered, or 'center_element' is reached.

        Notes
        -----
        Supports plane stress, plane strain, and 3D continua. Reduced plane
        strain requires sigma_zz = nu*(sigma_xx + sigma_yy), restricted to
        the built-in small-strain isotropic elastic law without initial
        stress or strain. Averaging responses differs from first averaging
        the stress tensors and then computing their invariant.

        Examples
        --------
        With a solved model part:

        >>> stresses = ResponseExtractor().extract_von_mises(model_part)
        """
        return self._extract_response(model_part, mode, self._von_mises)

    def extract_signed_von_mises(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> np.ndarray:
        """Extract von Mises stress with a tension/compression reporting sign.

        Mathematical Formulation
        ------------------------
        j = argmax_i(abs(sigma_i))
        signed_sigma_vm = sign(sigma_j) * sqrt(3/2 * sum_ij(s_ij**2))

        The sigma_i are descending principal Cauchy stresses and
        s = sigma - trace(sigma)/3 * I. Opposite-sign principal magnitude
        ties favor the positive principal stress; zero magnitude stays zero.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized continuum elements and a computed
            solution state, subject to the limits of ``extract_von_mises``.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': arithmetic mean of signed GP responses,
            allowing cancellation. 'Gp_max' selects the largest response
            magnitude and retains its sign; ties select the first GP.
            'center_gp' selects the mapped reference-element centroid, or
            warns and calls the unimplemented 'center_element' fallback.

        Returns
        -------
        numpy.ndarray, shape (n_elements,)
            Signed stresses in model stress units and element iteration order.
            An empty model part returns an empty array.

        Raises
        ------
        ValueError
            If the mode or a stress-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support or recovery check fails, or 'center_element' is
            reached directly or through the no-center fallback.
        numpy.linalg.LinAlgError
            If the symmetric eigenvalue calculation does not converge.

        Notes
        -----
        The sign follows the dominant principal stress, not the stress trace.
        Stress formats and plane-strain recovery match ``extract_von_mises``.
        """
        return self._extract_response(
            model_part, mode, self._signed_von_mises, gp_max_absolute=True
        )

    def extract_tresca(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> np.ndarray:
        """Extract one Tresca equivalent stress per continuum element.

        Mathematical Formulation
        ------------------------
        sigma_T = sigma1 - sigma3 = 2*tau_max

        Principal Cauchy stresses satisfy sigma1 >= sigma2 >= sigma3.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized continuum elements and a computed
            solution state, subject to the limits of ``extract_von_mises``.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': arithmetic mean of GP Tresca values.
            'Gp_max' takes the largest GP value. 'center_gp' selects the
            mapped reference-element centroid, or warns and calls the
            unimplemented 'center_element' fallback.

        Returns
        -------
        numpy.ndarray, shape (n_elements,)
            Nonnegative Tresca stresses in model stress units and element
            iteration order. An empty model part returns an empty array.

        Raises
        ------
        ValueError
            If the mode or a stress-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support or recovery check fails, or 'center_element' is
            reached directly or through the no-center fallback.
        numpy.linalg.LinAlgError
            If the symmetric eigenvalue calculation does not converge.

        Notes
        -----
        Stress formats and plane-strain recovery match ``extract_von_mises``.
        """
        return self._extract_response(model_part, mode, self._tresca)

    def extract_maximum_shear_stress(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> np.ndarray:
        """Extract the maximum shear stress magnitude per continuum element.

        Mathematical Formulation
        ------------------------
        tau_max = (sigma1 - sigma3)/2

        Principal Cauchy stresses satisfy sigma1 >= sigma2 >= sigma3.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized continuum elements and a computed
            solution state, subject to the limits of ``extract_von_mises``.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': arithmetic mean of GP maximum-shear values.
            'Gp_max' takes the largest GP value. 'center_gp' selects the
            mapped reference-element centroid, or warns and calls the
            unimplemented 'center_element' fallback.

        Returns
        -------
        numpy.ndarray, shape (n_elements,)
            Nonnegative maximum shear stresses in model stress units and
            element iteration order. Empty model parts return empty arrays.

        Raises
        ------
        ValueError
            If the mode or a stress-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support or recovery check fails, or 'center_element' is
            reached directly or through the no-center fallback.
        numpy.linalg.LinAlgError
            If the symmetric eigenvalue calculation does not converge.

        Notes
        -----
        Stress formats and plane-strain recovery match ``extract_von_mises``.
        """
        return self._extract_response(
            model_part, mode, self._maximum_shear_stress
        )

    def extract_signed_maximum_shear_stress(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> np.ndarray:
        """Extract maximum shear with a tension/compression reporting sign.

        Mathematical Formulation
        ------------------------
        j = argmax_i(abs(sigma_i))
        signed_tau_max = sign(sigma_j) * (sigma1 - sigma3)/2

        Principal Cauchy stresses satisfy sigma1 >= sigma2 >= sigma3.
        Opposite-sign principal magnitude ties favor the positive stress.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized continuum elements and a computed
            solution state, subject to the limits of ``extract_von_mises``.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': arithmetic mean of signed GP values,
            allowing cancellation. 'Gp_max' selects the largest magnitude
            and retains its sign; ties select the first GP. 'center_gp'
            selects the mapped reference-element centroid, or warns and
            calls the unimplemented 'center_element' fallback.

        Returns
        -------
        numpy.ndarray, shape (n_elements,)
            Signed maximum shear stresses in model stress units and element
            iteration order. An empty model part returns an empty array.

        Raises
        ------
        ValueError
            If the mode or a stress-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support or recovery check fails, or 'center_element' is
            reached directly or through the no-center fallback.
        numpy.linalg.LinAlgError
            If the symmetric eigenvalue calculation does not converge.

        Notes
        -----
        This is a reporting sign, not a traction sign on a specified plane.
        Stress formats and plane-strain recovery match ``extract_von_mises``.
        """
        return self._extract_response(
            model_part, mode, self._signed_maximum_shear_stress,
            gp_max_absolute=True
        )

    def extract_principle_stresses(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> List[List[float]]:
        """Extract descending principal Cauchy stresses per continuum element.

        Mathematical Formulation
        ------------------------
        det(sigma - lambda*I) = 0
        sigma1 >= sigma2 >= sigma3

        The three eigenvalues include compressive and out-of-plane stress.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized continuum elements and a computed
            solution state, subject to the limits of ``extract_von_mises``.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': average each ordered principal value by
            rank across GPs. 'Gp_max' returns the complete triplet from the
            GP with the largest sigma1; ties select the first GP. Components
            are never maximized independently. 'center_gp' selects the mapped
            reference-element centroid, or warns and calls the unimplemented
            'center_element' fallback.

        Returns
        -------
        list of list of float
            One [sigma1, sigma2, sigma3] list per element in model stress
            units. The outer list follows element iteration order and is
            empty for an empty model part. Directions are not returned.

        Raises
        ------
        ValueError
            If the mode or a stress-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support or recovery check fails, or 'center_element' is
            reached directly or through the no-center fallback.
        numpy.linalg.LinAlgError
            If the symmetric eigenvalue calculation does not converge.

        Notes
        -----
        The public spelling 'principle' is retained for compatibility.
        Stress formats and plane-strain recovery match ``extract_von_mises``.
        """
        return self._extract_response(
            model_part, mode, self._principal_stresses
        ).tolist()

    def extract_principle_strains(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> List[List[float]]:
        """Extract descending principal values of the supported small strain.

        Mathematical Formulation
        ------------------------
        det(epsilon - lambda*I) = 0
        epsilon1 >= epsilon2 >= epsilon3

        Voigt engineering shear is converted as epsilon_ij = gamma_ij/2;
        supplied full tensors already contain tensor shear.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized, supported small-displacement
            elements and a computed solution state. Strain is dimensionless.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': average each ordered principal value by
            rank across GPs. 'Gp_max' returns the complete triplet from the
            GP with the largest epsilon1; ties select the first GP.
            'center_gp' selects the mapped reference-element centroid, or
            warns and calls the unimplemented 'center_element' fallback.

        Returns
        -------
        list of list of float
            Dimensionless [epsilon1, epsilon2, epsilon3] lists in element
            iteration order. Empty model parts return an empty list.
            Principal directions are not returned.

        Raises
        ------
        ValueError
            If the mode or a strain-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support check fails, or 'center_element' is reached
            directly or through the no-center fallback.
        numpy.linalg.LinAlgError
            If the symmetric eigenvalue calculation does not converge.

        Notes
        -----
        Reads GREEN_LAGRANGE_STRAIN_VECTOR, evaluated with infinitesimal
        kinematics by the supported SmallDisplacementElement and 3D Bbar
        formulations. Only the exact built-in LinearElastic3DLaw,
        LinearElasticPlaneStress2DLaw, and LinearElasticPlaneStrain2DLaw
        types are accepted; finite-strain formulations are excluded.
        Missing epsilon_zz is zero for reduced plane strain by kinematics
        and for reduced plane stress by the Abaqus/CAE stored-field
        convention, even when physical thickness strain is nonzero.
        Explicit out-of-plane components are retained. The public spelling
        'principle' is retained for compatibility.
        """
        return self._extract_response(
            model_part, mode, self._principal_strains, strain_response=True
        ).tolist()

    def extract_von_mises_strain(
        self, model_part: KM.ModelPart, mode: str = "average"
    ) -> np.ndarray:
        """Extract the implemented Mises tensor invariant of small strain.

        Mathematical Formulation
        ------------------------
        epsilon_dev = epsilon - trace(epsilon)/3 * I
        epsilon_M = sqrt(3/2 * sum_ij(epsilon_dev_ij**2))

        Engineering shear is converted as epsilon_ij = gamma_ij/2 before
        evaluating the invariant. Full tensors already contain tensor shear.

        Parameters
        ----------
        model_part : KratosMultiphysics.ModelPart
            Model part with initialized, supported small-displacement
            elements and a computed solution state. Strain is dimensionless.
        mode : {'average', 'Gp_max', 'center_gp', 'center_element'}, optional
            Default is 'average': arithmetic mean of GP invariants.
            'Gp_max' takes the largest GP value. 'center_gp' selects the
            mapped reference-element centroid, or warns and calls the
            unimplemented 'center_element' fallback.

        Returns
        -------
        numpy.ndarray, shape (n_elements,)
            Nonnegative, dimensionless invariants in element iteration order.
            An empty model part returns an empty array.

        Raises
        ------
        ValueError
            If the mode or a strain-array shape is unsupported.
        RuntimeError
            If GP output is empty or its center-mode count is inconsistent.
        NotImplementedError
            If a support check fails, or 'center_element' is reached
            directly or through the no-center fallback.

        Notes
        -----
        Uses the same small-displacement linear-elastic support and reduced
        zero-epsilon_zz field convention as ``extract_principle_strains``.
        This applies the generic Abaqus/CAE MISES tensor operation to strain;
        it has not been numerically cross-validated against Abaqus. It is
        neither accumulated plastic strain PEEQ nor the alternative
        sqrt(2/3 * epsilon_dev:epsilon_dev) equivalent-strain definition.

        References
        ----------
        Abaqus FieldOutput, MISES and engineering-tensor handling:
        https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKERRefMap/simaker-c-fieldoutputpyc.htm
        """
        return self._extract_response(
            model_part, mode, self._von_mises_strain, strain_response=True
        )

    def _extract_response(
        self, model_part: KM.ModelPart, mode: str, response_measure,
        gp_max_absolute: bool = False, strain_response: bool = False
    ) -> np.ndarray:
        """Evaluate a response at GPs, then reduce it once per element.

        ``response_measure`` receives GP data, the analysis type, and (for
        stress) the optional recovery Poisson ratio. It returns a scalar or
        a descending principal-value triplet. ``gp_max_absolute`` selects
        signed scalars by magnitude; ``strain_response`` selects strain data
        and its stricter support checks. Both flags default to False.

        The result is an array of shape (n_elements,) or (n_elements, 3),
        except that all empty results have shape (0,). Validate and request
        data in the existing order, including before center-mode failures.
        """
        if mode not in ("center_gp", "average", "center_element", "Gp_max"):
            raise ValueError(
                "mode must be 'center_gp', 'average', 'center_element', "
                "or 'Gp_max'."
            )

        element_values = []
        for element in model_part.Elements:
            # Check supported kinematics before selecting or requesting data.
            analysis_type = self._get_analysis_type(element)
            if strain_response:
                self._check_strain_support(element)

            if mode == "center_element":
                element_values.append(self._extract_center_element(element))
                continue

            if mode == "center_gp":
                center_gp_index = self._get_center_gp_index(element)
                if center_gp_index is None:
                    KM.Logger.PrintWarning(
                        "ResponseExtractor",
                        f"Element {element.Id} has no Gauss point "
                        "at its center; "
                        "falling back to 'center_element'."
                    )
                    element_values.append(
                        self._extract_center_element(element)
                    )
                    continue

            # A plane-strain vector may retain sigma_zz that a 2D tensor omits.
            if strain_response:
                response_variable = KM.GREEN_LAGRANGE_STRAIN_VECTOR
            elif analysis_type == "plane_strain":
                response_variable = KM.CAUCHY_STRESS_VECTOR
            else:
                response_variable = KM.CAUCHY_STRESS_TENSOR

            gp_data = element.CalculateOnIntegrationPoints(
                response_variable, model_part.ProcessInfo
            )
            if not gp_data:
                raise RuntimeError(
                    f"Element {element.Id} returned no Gauss-point "
                    f"{'strain' if strain_response else 'stress'} tensors."
                )
            if (
                mode == "center_gp"
                and len(gp_data) != len(element.GetIntegrationPoints())
            ):
                raise RuntimeError(
                    f"Element {element.Id}: Gauss-point response count does "
                    "not match the active integration-point count."
                )

            # Recover omitted plane-strain stress only under the guarded law.
            poisson_ratio = None
            reduced_output = any(
                np.asarray(data).shape in ((3,), (2, 2)) for data in gp_data
            )
            if (
                reduced_output
                and not strain_response
                and analysis_type == "plane_strain"
            ):
                poisson_ratio = self._plane_strain_poisson_ratio(
                    element, model_part.ProcessInfo
                )

            if strain_response:
                response_arguments = (analysis_type,)
            else:
                response_arguments = (analysis_type, poisson_ratio)

            # Evaluate responses before averaging or selecting a GP maximum.
            if mode in ("average", "Gp_max"):
                gp_responses = []
                for data in gp_data:
                    gp_response = response_measure(data, *response_arguments)
                    gp_responses.append(gp_response)
                gp_responses = np.asarray(gp_responses)

                if mode == "average":
                    # r_element = (1/n_gp) * sum_g r_g, without GP weights.
                    # For principal values, average each already sorted rank.
                    element_response = np.mean(gp_responses, axis=0)
                elif gp_responses.ndim == 2:
                    # g* = argmax_g r_g[0]; retain the complete GP triplet.
                    element_response = gp_responses[
                        np.argmax(gp_responses[:, 0])
                    ]
                elif gp_max_absolute:
                    # g* = argmax_g |r_g|; keep its sign and first-GP ties.
                    element_response = gp_responses[
                        np.argmax(np.abs(gp_responses))
                    ]
                else:
                    element_response = np.max(gp_responses)
            else:
                element_response = response_measure(
                    gp_data[center_gp_index], *response_arguments
                )
            element_values.append(element_response)

        return np.asarray(element_values)

    @staticmethod
    def _check_strain_support(element: KM.Element) -> None:
        """Require supported small-displacement formulations and exact laws.

        Return None on success; otherwise raise NotImplementedError.
        Exact type checks intentionally exclude derived constitutive laws.
        """
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(
            element
        )
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        supported_formulation = (
            name.startswith("SmallDisplacementElement")
            or name.startswith("SmallDisplacementBbarElement3D")
        )
        if (
            supported_formulation
            and type(law) in (
                SMA.LinearElastic3DLaw,
                SMA.LinearElasticPlaneStress2DLaw,
                SMA.LinearElasticPlaneStrain2DLaw,
            )
        ):
            return
        raise NotImplementedError(
            f"Element {element.Id} ({name}): strain extraction currently "
            "requires a small-displacement continuum with a built-in "
            "linear-elastic law."
        )

    @staticmethod
    def _get_analysis_type(element: KM.Element) -> str:
        """Return 'plane_stress', 'plane_strain', or 'solid' from law features.

        Require matching local and working dimensions (2 or 3). Unsupported
        geometries or laws raise NotImplementedError. Plane-stress flags
        take precedence if a 2D law declares both plane flags.
        """
        geometry = element.GetGeometry()
        dimension = geometry.LocalSpaceDimension()
        if (
            dimension in (2, 3)
            and dimension == geometry.WorkingSpaceDimension()
            and element.Properties.Has(KM.CONSTITUTIVE_LAW)
        ):
            law = element.Properties[KM.CONSTITUTIVE_LAW]
            features = KM.ConstitutiveLawFeatures()
            law.GetLawFeatures(features)
            options = features.GetOptions()
            if dimension == 2:
                if options.Is(KM.ConstitutiveLaw.PLANE_STRESS_LAW):
                    return "plane_stress"
                if options.Is(KM.ConstitutiveLaw.PLANE_STRAIN_LAW):
                    return "plane_strain"
            elif (
                law.WorkingSpaceDimension() == 3
                and law.GetStrainSize() == 6
            ):
                return "solid"
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(
            element
        )
        raise NotImplementedError(
            f"Element {element.Id} ({name}): unsupported element or analysis; "
            "only plane stress, plane strain and 3D solids are supported."
        )

    @staticmethod
    def _plane_strain_poisson_ratio(
        element: KM.Element, process_info: KM.ProcessInfo
    ) -> float:
        """Return dimensionless nu for guarded isotropic stress recovery.

        The relation sigma_zz = nu*(sigma_xx + sigma_yy) requires a supported
        small-displacement element, exactly LinearElasticPlaneStrain2DLaw,
        and zero initial stress and strain at every GP. Missing initial-state
        output also rejects recovery. Raise NotImplementedError otherwise.
        """
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(
            element
        )
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        if (
            name.startswith((
                "SmallDisplacementElement", "SmallDisplacementBbarElement"
            ))
            and type(law) is SMA.LinearElasticPlaneStrain2DLaw
        ):
            for variable in (
                KM.INITIAL_STRESS_VECTOR, KM.INITIAL_STRAIN_VECTOR
            ):
                values = element.CalculateOnIntegrationPoints(
                    variable, process_info
                )
                if not values or any(
                    np.any(np.asarray(value) != 0.0) for value in values
                ):
                    break
            else:
                # Both initial-state checks passed without a break.
                return element.Properties[KM.POISSON_RATIO]
        raise NotImplementedError(
            f"Element {element.Id} ({name}): plane-strain output is missing "
            "sigma_zz. Recovery requires the built-in small-strain isotropic "
            "elastic law without initial stress or strain."
        )

    @staticmethod
    def _get_center_gp_index(element: KM.Element) -> Optional[int]:
        """Return the catalogued zero-based center index, or None if absent.

        Match the registered element and geometry, then the active GP count
        within the standard GI_GAUSS_1..5 family. The center is the mapped
        reference-element centroid, not necessarily Geometry.Center().
        Unknown elements or counts raise NotImplementedError. Index 0 is
        valid. Custom rules sharing a standard count cannot be distinguished.
        """
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(
            element
        )
        geometry_type = (
            element.GetGeometry().GetGeometryType().name
            .removeprefix("Kratos_")
        )
        element_data = _CENTER_GAUSS_POINT_LOOKUP["elements"].get(
            geometry_type
        )
        if element_data is None or name not in element_data["element_names"]:
            raise NotImplementedError(
                f"Element {element.Id} ({name}, {geometry_type}): no "
                "center-GP lookup for this element. Only the listed standard "
                "geometry-based GI_GAUSS_1 through GI_GAUSS_5 formulations "
                "are supported."
            )
        number_of_points = len(element.GetIntegrationPoints())
        rules = _CENTER_GAUSS_POINT_LOOKUP["families"][
            element_data["family"]
        ]["rules_by_gp_count"]
        rule = rules.get(number_of_points)
        if rule is None:
            raise NotImplementedError(
                f"Element {element.Id} ({name}): {number_of_points} "
                "integration points do not match a catalogued GI_GAUSS rule "
                "for this geometry."
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

    def _extract_center_element(self, element: KM.Element) -> NoReturn:
        """Reserve center recovery; always raise NotImplementedError."""
        raise NotImplementedError("'center_element' is not implemented yet.")

    @staticmethod
    def _von_mises(
        stress, analysis_type: str = "solid",
        poisson_ratio: Optional[float] = None
    ) -> np.floating:
        """Compute sqrt(3/2 * s:s) in the input stress units.

        Convert tensor or Kratos Voigt input using the analysis type and
        optional dimensionless recovery Poisson ratio. Here
        s = sigma - trace(sigma)/3 * I and s:s = sum_ij(s_ij**2).

        See Abaqus Conventions, Stress Invariants:
        https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMODRefMap/simamod-c-conventions.htm
        """
        tensor = ResponseExtractor._stress_tensor(
            stress, analysis_type, poisson_ratio
        )
        # Deviatoric Cauchy stress: s = sigma - trace(sigma)/3 * I.
        deviatoric_stress = tensor - np.trace(tensor) / 3.0 * np.eye(3)
        # Von Mises stress: sigma_vm = sqrt(3/2 * sum_ij(s_ij**2)).
        return np.sqrt(1.5 * np.sum(deviatoric_stress**2))

    @staticmethod
    def _signed_von_mises(
        stress, analysis_type: str = "solid",
        poisson_ratio: Optional[float] = None
    ) -> np.floating:
        """Return sign(sigma_j)*sigma_vm, with j = argmax_i |sigma_i|.

        Units and conversion arguments match ``_von_mises``. Descending
        principal ordering makes opposite-sign magnitude ties favor tension.
        """
        principal_stresses = ResponseExtractor._principal_stresses(
            stress, analysis_type, poisson_ratio
        )
        # j = argmax_i |sigma_i|; argmax retains the first ordered tie.
        dominant_stress = principal_stresses[
            np.argmax(np.abs(principal_stresses))
        ]
        # Signed invariant: signed_sigma_vm = sign(sigma_j) * sigma_vm.
        return np.sign(dominant_stress) * ResponseExtractor._von_mises(
            stress, analysis_type, poisson_ratio
        )

    @staticmethod
    def _tresca(
        stress, analysis_type: str = "solid",
        poisson_ratio: Optional[float] = None
    ) -> np.floating:
        """Return sigma1 - sigma3 in input stress units.

        Conversion arguments match ``_von_mises``. The symmetric eigensolver
        returns ascending values here, so the last minus the first is used.
        """
        tensor = ResponseExtractor._stress_tensor(
            stress, analysis_type, poisson_ratio
        )
        principal_stresses = np.linalg.eigvalsh(tensor)
        # Tresca equivalent stress: sigma_T = sigma_max - sigma_min.
        return principal_stresses[-1] - principal_stresses[0]

    @staticmethod
    def _maximum_shear_stress(
        stress, analysis_type: str = "solid",
        poisson_ratio: Optional[float] = None
    ) -> np.floating:
        """Return (sigma1 - sigma3)/2 using ``_tresca`` arguments and units."""
        # Maximum shear stress: tau_max = sigma_T/2.
        return 0.5 * ResponseExtractor._tresca(
            stress, analysis_type, poisson_ratio
        )

    @staticmethod
    def _signed_maximum_shear_stress(
        stress, analysis_type: str = "solid",
        poisson_ratio: Optional[float] = None
    ) -> np.floating:
        """Return sign(sigma_j)*(sigma1 - sigma3)/2 in input stress units.

        Conversion arguments match ``_von_mises``. The index j maximizes
        |sigma_j|; opposite-sign principal magnitude ties favor tension.
        """
        principal_stresses = ResponseExtractor._principal_stresses(
            stress, analysis_type, poisson_ratio
        )
        # j = argmax_i |sigma_i|; descending order favors a positive tie.
        dominant_stress = principal_stresses[
            np.argmax(np.abs(principal_stresses))
        ]
        # Signed maximum shear: sign(sigma_j) * (sigma1 - sigma3)/2.
        return np.sign(dominant_stress) * 0.5 * (
            principal_stresses[0] - principal_stresses[-1]
        )

    @staticmethod
    def _principal_stresses(
        stress, analysis_type: str = "solid",
        poisson_ratio: Optional[float] = None
    ) -> np.ndarray:
        """Return a descending three-value stress array in input stress units.

        Conversion arguments match ``_von_mises``. Symmetric tensor input is
        expected; no symmetrization or eigenvector calculation is performed.
        """
        tensor = ResponseExtractor._stress_tensor(
            stress, analysis_type, poisson_ratio
        )
        # det(sigma - lambda*I) = 0; reverse ascending symmetric eigenvalues.
        return np.linalg.eigvalsh(tensor)[::-1]

    @staticmethod
    def _principal_strains(
        strain, analysis_type: str = "solid"
    ) -> np.ndarray:
        """Return a descending three-value array of dimensionless strains.

        Convert the input with ``_strain_tensor``. Symmetric tensors are
        expected; no symmetrization or eigenvector calculation is performed.
        """
        tensor = ResponseExtractor._strain_tensor(strain, analysis_type)
        # det(epsilon - lambda*I) = 0; epsilon1 >= epsilon2 >= epsilon3.
        return np.linalg.eigvalsh(tensor)[::-1]

    @staticmethod
    def _von_mises_strain(
        strain, analysis_type: str = "solid"
    ) -> np.floating:
        """Return sqrt(3/2 * epsilon_dev:epsilon_dev), dimensionless.

        Input conversion follows ``_strain_tensor``. This is the implemented
        Mises tensor invariant, not accumulated equivalent plastic strain.
        """
        tensor = ResponseExtractor._strain_tensor(strain, analysis_type)
        # Deviatoric strain: epsilon_dev = epsilon - trace(epsilon)/3 * I.
        deviatoric_strain = tensor - np.trace(tensor) / 3.0 * np.eye(3)
        # Mises strain invariant: sqrt(3/2 * sum_ij(epsilon_dev_ij**2)).
        return np.sqrt(1.5 * np.sum(deviatoric_strain**2))

    @staticmethod
    def _strain_tensor(strain, analysis_type: str = "solid") -> np.ndarray:
        """Convert dimensionless tensor/Voigt strain to a 3x3 tensor.

        Full tensors retain supplied components. Six-component Voigt order
        is [xx, yy, zz, gamma_xy, gamma_yz, gamma_xz]. Plane analyses also
        accept [xx, yy, gamma_xy], [xx, yy, zz, gamma_xy], and 2x2 tensors.
        Divide Voigt engineering shear by two, but never halve tensor shear.
        Missing epsilon_zz is zero by plane-strain kinematics or the reduced
        plane-stress stored-field convention. Full 3x3 arrays may be returned
        without a copy. Unsupported analyses raise NotImplementedError;
        incompatible shapes raise ValueError.

        https://docs.software.vt.edu/abaqusv2025/English/SIMACAEOUTRefMap/simaout-c-stdoutputvar.htm
        """
        strain = np.asarray(strain)
        if analysis_type not in ("plane_stress", "plane_strain", "solid"):
            raise NotImplementedError(
                f"Unsupported analysis type: {analysis_type}."
            )
        if strain.shape == (3, 3):
            return strain
        if strain.shape == (6,):
            xx, yy, zz, gamma_xy, gamma_yz, gamma_xz = strain
            # Tensor shear strain: epsilon_ij = gamma_ij/2 for i != j.
            return np.array([
                [xx, gamma_xy / 2.0, gamma_xz / 2.0],
                [gamma_xy / 2.0, yy, gamma_yz / 2.0],
                [gamma_xz / 2.0, gamma_yz / 2.0, zz],
            ])
        if analysis_type in ("plane_stress", "plane_strain"):
            if strain.shape == (4,):
                xx, yy, zz, gamma_xy = strain
                # Preserve supplied epsilon_zz; epsilon_xy = gamma_xy/2.
                return np.array([
                    [xx, gamma_xy / 2.0, 0.0],
                    [gamma_xy / 2.0, yy, 0.0],
                    [0.0, 0.0, zz],
                ])
            tensor = np.zeros((3, 3))
            if strain.shape == (2, 2):
                tensor[:2, :2] = strain
            elif strain.shape == (3,):
                xx, yy, gamma_xy = strain
                # Only Voigt input stores engineering shear gamma_xy.
                tensor[:2, :2] = [
                    [xx, gamma_xy / 2.0],
                    [gamma_xy / 2.0, yy],
                ]
            else:
                raise ValueError(
                    f"Unsupported {analysis_type} strain shape: "
                    f"{strain.shape}."
                )
            # epsilon_zz = 0 for reduced plane strain by kinematics and for
            # reduced plane stress by the stored-field convention; no physical
            # thickness-strain recovery is performed.
            return tensor
        raise ValueError(
            "A 3D solid requires a 3x3 strain tensor "
            "or six-component Voigt vector."
        )

    @staticmethod
    def _stress_tensor(
        stress, analysis_type: str, poisson_ratio: Optional[float]
    ) -> np.ndarray:
        """Convert tensor/Voigt stress to a 3x3 tensor in input stress units.

        Six-component Kratos order is [xx, yy, zz, xy, yz, xz], with physical
        shear stresses. Plane analyses also accept [xx, yy, xy],
        [xx, yy, zz, xy], and 2x2 tensors. Preserve full supplied components;
        missing plane-stress components are zero. Reduced plane strain needs
        dimensionless poisson_ratio for sigma_zz = nu*(sigma_xx + sigma_yy),
        guarded by the caller. Full 3x3 arrays may be returned without a copy.
        Unsupported analyses or missing recovery data raise
        NotImplementedError; incompatible shapes raise ValueError.
        """
        stress = np.asarray(stress)
        if analysis_type not in ("plane_stress", "plane_strain", "solid"):
            raise NotImplementedError(
                f"Unsupported analysis type: {analysis_type}."
            )
        if stress.shape == (3, 3):
            tensor = stress
        elif stress.shape == (6,):
            # Kratos order: xx, yy, zz, xy, yz, xz; no shear scaling.
            xx, yy, zz, xy, yz, xz = stress
            tensor = np.array([
                [xx, xy, xz],
                [xy, yy, yz],
                [xz, yz, zz],
            ])
        elif analysis_type in ("plane_stress", "plane_strain"):
            tensor = np.zeros((3, 3))
            if stress.shape == (4,):
                # Kratos ordering: xx, yy, zz, xy.
                xx, yy, zz, xy = stress
                tensor = np.array([
                    [xx, xy, 0.0],
                    [xy, yy, 0.0],
                    [0.0, 0.0, zz],
                ])
            else:
                if stress.shape == (2, 2):
                    tensor[:2, :2] = stress
                elif stress.shape == (3,):
                    xx, yy, xy = stress
                    tensor[:2, :2] = [
                        [xx, xy],
                        [xy, yy],
                    ]
                else:
                    raise ValueError(
                        f"Unsupported {analysis_type} stress shape: "
                        f"{stress.shape}."
                    )
                # Omitted plane-stress components are zero by definition.
                # Plane strain requires sigma_zz from the material response.
                if analysis_type == "plane_strain":
                    if poisson_ratio is None:
                        raise NotImplementedError(
                            "Plane-strain stress extraction requires sigma_zz."
                        )
                    # Isotropic elastic recovery: sigma_zz = nu*(sigma_xx
                    # + sigma_yy). The caller checks the law and initial state.
                    tensor[2, 2] = poisson_ratio * (
                        tensor[0, 0] + tensor[1, 1]
                    )
        else:
            raise ValueError(
                "A 3D solid requires a 3x3 tensor "
                "or six-component Voigt vector."
            )
        return tensor
