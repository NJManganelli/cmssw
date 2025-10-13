import argparse
import correctionlib
import correctionlib.schemav2 as cs
import json
import rich
import copy

def build_relative_and_compound_corrections(in_file, out_file, override_flow="clamp", toysim_2d_activate_z0_swizzle=False):
    """
    Take SmartPixels ToyDetector pt-parameterized resolutions and build relative corrections
    between configuration and the no-SmartPixels ('0000' / 'IIII') resolutions.
    Then, build compound corrections that combine the absolute resolutions with HashPRNG nodes
    and generalize inputs to include pt, eta, phi, z0, and d0 of tracking particle (TP). Outputs absolute smearing to add to TP.
    Then build compound corrections combining the relative resolutions with actual reconstructed versus TP kinematics,
    generalizing to TP kinematics, TP-to-track matching category, and the reconstructed L1Track kinematics. 
    Outputs absolute smearing to add to TP.

    toysim_2d_activate_z0_swizzle should be utilized when 2D toysim files are converted, in order to create a z0 resolution copying the d0 relative resolution but the z0 absolute difference
    """
    to_replace = ["_0001", "_0010", "_0100", "_1000",
                  "_0011", "_0101", "_1001", "_0110", "_1010", "_1100",
                  "_0111", "_1011", "_1101", "_1110",
                  "_1111"]
    cset = correctionlib.CorrectionSet.from_file(in_file)
    json_dict = None
    with open(in_file, 'rb') as raw_data:
        json_dict = json.load(raw_data)
        json_dict = {cdict['name']: cdict for cdict in json_dict['corrections']}


    hash_prng_node = cs.Correction(
        name="hash_prng_node",
        description="Node to generate a random seed based on tracking particle kinematics",
        version=1,
        inputs = [
            cs.Variable(name="pt_tp", type="real", description="Tracking particle transverse momentum"),
            cs.Variable(name="eta_tp", type="real", description="Tracking particle pseudorapidity"),
            cs.Variable(name="phi_tp", type="real", description="Tracking particle azimuthal angle"),
            cs.Variable(name="z0_tp", type="real", description="Tracking particle signed longitudinal impact parameter"),
            cs.Variable(name="d0_tp", type="real", description="Tracking particle signed transverse impact parameter"),
        ],
        output=cs.Variable(name="rng", type="real"),
        data=cs.HashPRNG(
            nodetype="hashprng",
            inputs=["pt_tp", "eta_tp", "phi_tp", "z0_tp", "d0_tp"],
            distribution="stdnormal",
        )
    )
    new_cset_corrections = [hash_prng_node]
    for variable in ["pt", "eta", "phi", "z0", "d0"]:
        new_cset_corrections.append(
            cs.Correction(
                name=f"{variable}_track_tp_difference",
                version=1,
                description=f"L1Track minus Tracking Particle {variable} value",
                inputs=[
                    cs.Variable(name=f"{variable}_track", type="real", description=f"L1 track {variable}"),
                    cs.Variable(name=f"{variable}_tp", type="real", description=f"Tracking particle {variable}"),
                ],
                output=cs.Variable(name=f"{variable}_track_tp_difference", type="real", description=f"L1 Track minus Tracking Particle {variable}"),
                data=cs.Formula(
                    nodetype="formula",
                    variables=[f"{variable}_track", f"{variable}_tp"], #[ x, y, z, t ]...
                    parser="TFormula",
                    expression="(x - y)",
                ),
            )
        )
    new_cset_compound_corrections = []

    if not any(["z0" in k for k in cset.keys()]):
        print("No z0 corrections found, pass toysim_2d_activate_z0_swizzle=True or --toysim_2d_activate_z0_swizzle to inser swizzled z0 compound corrections")
    for key in cset.keys():
        variable = key.split("_")[0]
        if variable not in ["pt", "eta", "phi", "z0", "d0"]:
            raise NotImplementedError(f"Variable {variable} is not supported for relative smears or compounding corrections")
        corr = cset[key]
        denominator_key = key
        for tag in to_replace:
            denominator_key = denominator_key.replace(tag, "_0000")
        corr_numerator_data = json_dict[key]['data']
        corr_denominator_data = json_dict[denominator_key]['data']

        if corr_numerator_data['nodetype'] != "binning":
            print(f"Skipping {key} as it is not a binning correction")
            continue
        relative_data = {}
        relative_data['nodetype'] = "binning"
        relative_data['input'] = corr_numerator_data['input']
        relative_data['edges'] = corr_numerator_data['edges']
        relative_data['content'] = []
        for i, corr_numerator_data_content_i in enumerate(corr_numerator_data['content']):
            if isinstance(corr_numerator_data_content_i, int):
                if corr_denominator_data['content'][i] == 0:
                    relative_data['content'].append(0)
                else:
                    relative_data['content'].append(corr_numerator_data_content_i / corr_denominator_data['content'][i])
            elif isinstance(corr_numerator_data_content_i, dict):
                corr_numerator_data_dict = copy.deepcopy(corr_numerator_data_content_i)
                corr_denominator_data_dict = corr_denominator_data['content'][i] # just a clear reference
                for j, corr_num_data_content_i_j in enumerate(corr_numerator_data_dict['content']):
                    if (corr_num_data_content_i_j == 0) and (corr_denominator_data_dict['content'][j] == 0):
                        # both are 0, this is our regime where coverage ends, e.g. |eta| > 0.8 in v2p0, so we return a ratio of 1 and permit clamping to it
                        corr_numerator_data_dict['content'][j] = 1.0
                    elif corr_denominator_data_dict['content'][j] == 0:
                        # just 0 denominator, though this code path should be inaccessible...
                        corr_numerator_data_dict['content'][j] = 0.0
                    else:
                        # Take the ratio
                        corr_numerator_data_dict['content'][j] = corr_num_data_content_i_j / corr_denominator_data_dict['content'][j]
                corr_numerator_data_dict['flow'] = override_flow if override_flow else corr_numerator_data_dict.get('flow')
                # now append to the (outer) list contents
                relative_data['content'].append(corr_numerator_data_dict)
            else:
                raise NotImplementedError(f"Unhandled level-1 content type: {corr_numerator_data_content_i}")
        relative_data['flow'] = override_flow if override_flow else corr_numerator_data.get('flow') #corr_numerator_data.get('flow', 'clamp')

        new_inputs = [cs.Variable(name=var.name, type=var.type, description=var.description) for var in corr.inputs]
        new_output = cs.Variable(name=corr.output.name, type=corr.output.type, description=corr.output.description)
        new_corr_key = key.replace("_smear", "_relative_smear")

        old_corr_reconstituted = cs.Correction(
            name=key,
            description=corr.description,
            version=corr.version,
            inputs=new_inputs,
            output=new_output,
            data=cs.Binning(
                nodetype=corr_numerator_data['nodetype'],
                input=corr_numerator_data['input'],
                edges=corr_numerator_data['edges'],
                content=corr_numerator_data['content'],
                flow=override_flow if override_flow else corr_numerator_data.get('flow')
            )
        )
        new_corr = cs.Correction(
            name=corr.name.replace("_smear", "_relative_smear"),
            description=corr.description,
            version=corr.version,
            #inputs=corr.inputs,
            #output=corr.output,
            #data = cs.Binning(**relative_data)
            inputs=new_inputs,
            output=new_output,
            data=cs.Binning(
                nodetype=relative_data['nodetype'],
                input=relative_data['input'],
                edges=relative_data['edges'],
                content=relative_data['content'],
                flow=relative_data['flow'],
            )
            #    nodetype="binning",
            #    input="pt",
            #    edges=[10, 20, 30, 40, 50, 80, 120],
            #    content=[0.3, 0.25, 0.20, 0.14, 0.06, 0.02],
            #    flow="clamp",
            #)
        )
        #new_corr.data.content[0].value.expression = "x + 0.1 * x"  # Example modification
        new_cset_corrections.append(old_corr_reconstituted) # copy the old correction too
        new_cset_corrections.append(new_corr) # add the new correction

        # Now create a compound correction using the old correction and a HashPRNG node
        new_cset_compound_corrections.append(
            cs.CompoundCorrection(
                name=key.replace("_smear", "_smear_compound"),
                description="",
                inputs=[
                    cs.Variable(name="pt_tp", type="real", description="Tracking particle transverse momentum"),
                    cs.Variable(name="eta_tp", type="real", description="Tracking particle pseudorapidity"),
                    cs.Variable(name="phi_tp", type="real", description="Tracking particle azimuthal angle"),
                    cs.Variable(name="z0_tp", type="real", description="Tracking particle signed longitudinal impact parameter"),
                    cs.Variable(name="d0_tp", type="real", description="Tracking particle signed transverse impact parameter"),
                ],
                output=cs.Variable(name="shift", type="real", description=f"Additive shift to Truth Particle {variable}"),
                inputs_update=[],
                input_op="*",
                output_op="*",
                stack=[key, "hash_prng_node"],
            )
        )

        # Now create the compound correction that combines the relative correction with a Formula shift node
        new_cset_compound_corrections.append(
            cs.CompoundCorrection(
                name=key.replace("_smear", "_relative_smear_compound"),
                description="",
                inputs=[
                    cs.Variable(name="pt_tp", type="real", description="Tracking particle transverse momentum"),
                    cs.Variable(name="eta_tp", type="real", description="Tracking particle pseudorapidity"),
                    cs.Variable(name="phi_tp", type="real", description="Tracking particle azimuthal angle"),
                    cs.Variable(name="z0_tp", type="real", description="Tracking particle signed longitudinal impact parameter"),
                    cs.Variable(name="d0_tp", type="real", description="Tracking particle signed transverse impact parameter"),
                    cs.Variable(name="tp_track_match", type="int", description="match class of the tracking particle to an L1 track (4: isLooselyGenuine, 12: isLooselyGenuine && isGenuine"),
                    cs.Variable(name="pt_track", type="real", description="Reconstructed L1 track transverse momentum"),
                    cs.Variable(name="eta_track", type="real", description="Reconstructed L1 track pseudorapidity"),
                    cs.Variable(name="phi_track", type="real", description="Reconstructed L1 track azimuthal angle"),
                    cs.Variable(name="z0_track", type="real", description="Reconstructed L1 track signed longitudinal impact parameter"),
                    cs.Variable(name="d0_track", type="real", description="Reconstructed L1 track signed transverse impact parameter"),
                ],
                output=cs.Variable(name="shift", type="real", description=f"Additive shift to Truth Particle {variable}"),
                inputs_update=[],
                input_op="*",
                output_op="*",
                stack=[new_corr_key, f"{variable}_track_tp_difference"],
            )
        )
        if "d0" in key and toysim_2d_activate_z0_swizzle:
            print("Creating z0 compound correction using d0 relative smear and z0_track_tp_difference for SPix config ", key.split("_")[-1])
            new_cset_compound_corrections.append(
                cs.CompoundCorrection(
                    name=key.replace("d0", "z0").replace("_smear", "_relative_smear_compound"),
                    description="",
                    inputs=[
                        cs.Variable(name="pt_tp", type="real", description="Tracking particle transverse momentum"),
                        cs.Variable(name="eta_tp", type="real", description="Tracking particle pseudorapidity"),
                        cs.Variable(name="phi_tp", type="real", description="Tracking particle azimuthal angle"),
                        cs.Variable(name="z0_tp", type="real", description="Tracking particle signed longitudinal impact parameter"),
                        cs.Variable(name="d0_tp", type="real", description="Tracking particle signed transverse impact parameter"),
                        cs.Variable(name="tp_track_match", type="int", description="match class of the tracking particle to an L1 track (-1: isFake, 0: isUnknown, 1: isLooselyGenuine, 2: isGenuine, 3: isCombinatoric)"),
                        cs.Variable(name="pt_track", type="real", description="Reconstructed L1 track transverse momentum"),
                        cs.Variable(name="eta_track", type="real", description="Reconstructed L1 track pseudorapidity"),
                        cs.Variable(name="phi_track", type="real", description="Reconstructed L1 track azimuthal angle"),
                        cs.Variable(name="z0_track", type="real", description="Reconstructed L1 track signed longitudinal impact parameter"),
                        cs.Variable(name="d0_track", type="real", description="Reconstructed L1 track signed transverse impact parameter"),
                    ],
                    output=cs.Variable(name="shift", type="real", description=f"Additive shift to Truth Particle {variable}"),
                    inputs_update=[],
                    input_op="*",
                    output_op="*",
                    stack=[new_corr_key, f"z0_track_tp_difference"],
                )
            )
        # if (key.endswith("1111") or key.endswith("0000")) and ("pt" in key or "d0" in key):
        #     rich.print(new_cset_compound_corrections[-2])
        #     rich.print(new_cset_compound_corrections[-1])

    new_cset = correctionlib.schemav2.CorrectionSet(
        schema_version=2,
        corrections=new_cset_corrections,
        compound_corrections=new_cset_compound_corrections,
    )
    with open(out_file, "w") as fout:
        fout.write(new_cset.model_dump_json(exclude_unset=True))

if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("--in_file" ,   type=str, default="spixel_smear_all_configs_labeled_json.json", help="input json file with absolute resolutions")
    parser.add_argument("--out_file" ,   type=str, default="spixel_smear_all_configs_labeled_json_compound_z0swizzle.json", help="output json file with absolute, relative, and compound resolutions")
    parser.add_argument("--override_flow", type=str, default="clamp", help="Override the input json flow behavior to this type, set to 'none' to not override")
    parser.add_argument("--toysim_2d_activate_z0_swizzle", action='store_true')
    options = parser.parse_args()

    override_flow = options.override_flow if options.override_flow.lower() != "none" else None
    build_relative_and_compound_corrections(options.in_file, options.out_file, override_flow=override_flow, toysim_2d_activate_z0_swizzle=options.toysim_2d_activate_z0_swizzle)
