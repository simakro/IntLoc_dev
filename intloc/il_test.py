# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
import sys
import json
from collections import defaultdict

from il_logging import Ilogger
from il_html import read_il_csv
from il_chk_conv import report_time
from il_subprocesses import intloc_self_sp


ilog = Ilogger()
ilog.module = __name__


def compare_sim_ints_to_results(
        species,
        outdir,
        int_rep_path,
        req_precision,
        report="Integration_Report.csv"
        ):
    """Compare intloc results to simulated integrations in test data"""
    report_time(
        "comparison of intloc results with record of simulated ints in testdata"
        )
    if any([species=="Saccharomyces cerevisiae", species=="yeast"]):
        integrator_report = int_rep_path
    elif any([species=="Homo sapiens", species=="human"]):
        integrator_report = int_rep_path
    else:
        ilog.vlprint(
            f"WARNING: No test data available for species {species}. Skipping"
            f" testing for correct positions of integrations.", 2
            )

    with open(integrator_report, "r") as sim_ints:
        # Here also extract ID of Integration
        correct_loc = defaultdict(list)
        check_list = []
        for line in sim_ints:
            try:
                ls = line.strip().split(",")
                id = int(ls[0])
                chrom = ls[1].strip().split(" ")[0]
                loc = int(ls[-1].strip())
                correct_loc[chrom].append((loc, id))
                check_list.append(chrom + "_" + str(loc))
            except ValueError:
                pass

    
    results_csv = os.path.join(outdir, report)
    with open(results_csv, "r") as il_results:
        # results_loc = dict()
        results_loc = defaultdict(list)
        for line in il_results:
            if line.startswith("ID"):
                pass
            else:
                ls = line.strip().split(",")
                chrom = ">" + ls[2].strip().split(" ")[0]
                loc = int(ls[3].strip())
                # if chrom in results_loc:
                results_loc[chrom].append(loc)
                # else:
                #     results_loc[chrom] = [loc]

    test_results = os.path.join(outdir, "test_results.csv")
    testres_dct = {}
    corr2found = {}
    redundant_sites = 0
    with open(test_results, "w") as tres:
        tres.write("Int_called, Nearest_correct, ID_nearest, diff[bp]\n")
        for chrom in results_loc:
            if chrom in correct_loc:
                for res_loc in results_loc[chrom]:
                    best_match = min(
                        [
                            (abs(res_loc-corr_loc[0]), corr_loc) for corr_loc in correct_loc[chrom]
                            ]
                            ,key=lambda x: x[0]
                        )
                    name_best_match = chrom + "_" + str(best_match[1][0])
                    if best_match[0] < req_precision:
                        # name_best_match = chrom + "_" + str(best_match[1][0])
                        try:
                            check_list.remove(name_best_match)
                            corr2found[name_best_match] = f"{chrom}_{res_loc}"
                        except ValueError:
                            ilog.vlprint(
                                f"Integration {name_best_match} has already"
                                f" been matched with another PotIntSite", 0
                            )
                    ilog.vlprint(
                        f"Result {chrom}_{res_loc} differs {best_match[0]} bp "
                        f"from best matching actual integration ID "
                        f"{best_match[1][1]} {name_best_match}", 0
                        # f"{best_match[1][1]} {chrom}_{best_match[1][0]}", 0
                        )
                    tres.write(
                        f"{chrom}_{res_loc},{name_best_match},"
                        f"{best_match[1][1]},{best_match[0]}\n"
                        )
                    if name_best_match not in testres_dct:
                        testres_dct[name_best_match] = best_match[0]
                    else:
                        redundant_sites += 1
                        if best_match[0] < testres_dct[name_best_match]:
                             testres_dct[name_best_match] = best_match[0]
            else:
                for res_loc in results_loc[chrom]:
                    print(
                        f"Result {chrom}_{res_loc} does not correspond to an " 
                        f"actual integration, since no integrations have been "
                        f"inserted into {chrom}"
                        )
    
        for int_name in check_list:
            ilog.vlprint(f"Integration {int_name} is missing from intloc results.", 0)
            tres.write(f"Integration {int_name} is missing from intloc results.\n")
            if int_name not in testres_dct:
                testres_dct[int_name] = None
    report_time("Testing of results from simulated data", end=True)
    return testres_dct, redundant_sites, corr2found


def report_int_comp_results(loc_res, redundant, tolerance):
    diffs_vals = loc_res.values()
    diffs = [d for d in diffs_vals if d!=None]
    diff_cats = {
        "perfect (0)": len([d for d in diffs if d==0]),
        "near_perf (1)": len([d for d in diffs if d==1]),
        "very_close (10)":  len([d for d in diffs if 1<d<11]),
        "close (25)":  len([d for d in diffs if 10<d<26]),
        "acceptable (100)":  len([d for d in diffs if 25<d<101]),
        "far (>100)":  len([d for d in diffs if 100<d<tolerance+1]),
        "out of tolerance": len([d for d in diffs if tolerance<d]),
        "avg. divergence": sum(diffs)/len(diffs),
        "not present": len([d for d in diffs_vals if d==None])
    }
    print("Matches summary:")
    for cat in diff_cats:
        if diff_cats[cat] > 0:
            print(f"\t{cat}: {diff_cats[cat]}")
    
    if diff_cats["not present"]==0 and diff_cats["out of tolerance"]==0:
        print(
            f"All int-sites found within req. precision of {tolerance}."
            )
    elif diff_cats["not present"]==0:
        print(
            f"All int-sites found, but {diff_cats['out of tolerance']}" 
            f" out of range {tolerance}."
            )
    else:
        print(f"{diff_cats['not present']} int-sites are missing.")

    if redundant>0:
        print(f"{redundant} reduntant/supernumerary sites were called")


def compare_run_info(report, param, expect, run_inf, kwargs):
    """Compare data in run_info.csv to expectations"""
    deviation_found = False
    print(f"Expected {param}: {expect[report][param]}")
    result = eval(run_inf[report][param])
    if result == expect[report][param]:
        print(f"Expected result was generated: {result}")
    else:
        print(f"Divergent result was generated {result}")
        deviation_found = True
    return deviation_found

                    # report, entry, expect, run_inf, kwargs
def compare_int_rep_cov(report, corr_site, expect, run_inf, kwargs):
    """Compare data in run_info.csv to expectations"""
    deviation_found = False
    # found[1:] to remove ">" from int-name
    sites = {k[1:]:v[1:] for k,v in kwargs["sites"].items()}
    int_rep = {s["INT_name"]:s for s in run_inf[report]}
    by_corr_site = {corr:int_rep[found] for corr,found in sites.items()}
    if corr_site in by_corr_site:
        print(f"Real int-site {corr_site}")
        for param in expect[report][corr_site]:
            xpct_val = expect[report][corr_site][param]
            print(f"\t{param} expected value {xpct_val}")
            found_val = int_rep[sites[corr_site]][param]
            try:
                found_val = int(found_val)
            except ValueError:
                pass
            if xpct_val == found_val:
                print(f"\tFound expected value.")
            else:
                print(f"\tFound deviating value {found_val}.")
                deviation_found = True
    else:
        print(f"No match was found for expected real int-site {corr_site}")
    return deviation_found


def run_batch_test(batch_file, outdir=False, tolerance=1000):
    """Run a batch of intloc commands and assess results"""
    with open(batch_file, "r") as data:
        runs = json.load(data)
        paths_keys = ["integrator_file", "integrator_file_preex"]
        # generate platform and install-dir specific paths for test-data
        for run in runs:
            for key in runs[run]:
                module_path = os.path.split(__file__)[0]
                if key in paths_keys and type(runs[run][key])==list:
                    runs[run][key].insert(0, module_path)
                    runs[run][key] = os.path.join(*runs[run][key])
            func_call = " ".join(runs[run]["cmd"]["pos_flag_args"])
            for kwarg in runs[run]["cmd"]["kw_args"]:
                func_call = " ".join([func_call, kwarg])
                value = runs[run]["cmd"]["kw_args"][kwarg]
                if type(value)==list:
                    value.insert(0, module_path)
                    value = os.path.join(*value)
                func_call = " ".join([func_call, str(value)])
            runs[run]["cmd"] = func_call
    
    for run_name in runs:
        output = intloc_self_sp(runs[run_name]["cmd"], outdir=outdir)
        out = output.stdout.decode().split("\n")
        out = [l for l in out[-5:] if l.startswith("Results saved to:")][0]
        out_dir = out.split("saved to:")[1].strip()
        runs[run_name]["out_dir"] = out_dir
        rep_dir = os.path.join(out_dir, "reports")
        for file in [f for f in os.scandir(rep_dir) if f.name.endswith(".csv")]:
            if file.name != "run_info.csv":
                inf = read_il_csv(file.path)
                runs[run_name][file.name] = inf
            else:
                with open(file.path, "r") as inf:
                    data = inf.read().split("\n")
                    data = [l for l in data if len(l.strip())>0]
                run_info = {k:v for k,v in [l.split(",") for l in data]}
                runs[run_name][file.name] = run_info
    # import pprint
    # pprint.pprint(runs[run]["seq_gain_loss_at_intsites.csv"])

    no_deviations_in_any_run = True
    for run in runs:
        print(f"##################{run}##########################")
        print("###############################################################")
        run_inf = runs[run]
        expect = run_inf["expectations"]
        mod_out = os.path.join(run_inf["out_dir"], "reports")

        loc_res, redundant, corr2act = compare_sim_ints_to_results(
            run_inf["species"], mod_out, run_inf["integrator_file"], tolerance
            )
        report_int_comp_results(loc_res, redundant, tolerance)
        if run_inf["mode"] == "intra":
            loc_res, redundant, c2a_intra_map = compare_sim_ints_to_results(
                run_inf["species"],
                mod_out,
                run_inf["integrator_file_preex"],
                tolerance,
                report="preex_Integration_Report.csv"
            )
            report_int_comp_results(loc_res, redundant, tolerance)
        
        comp_func = {
            "run_info.csv": (compare_run_info,{}),
            "Integration_Report.csv": (compare_int_rep_cov,{"sites": corr2act}),
        }
        deviation_found = False
        for report in expect:
            print(f"Comparing {report} to expectations")
            for entry in expect[report]:
                func, kwargs = comp_func[report]
                dev_found = func(report, entry, expect, run_inf, kwargs)
                if dev_found:
                    deviation_found = True
                # print(f"Expected {param}: {expect[report][param]}")
                # result = eval(run_inf[report][param])
                # if result == expect[report][param]:
                #     print(f"Expected result was generated: {result}")
                # else:
                #     print(f"Divergent result was generated {result}")
        
        if not deviation_found:
            print(f"No deviation from expectations was found in run {run}")
        else:
            print(f"Some deviation(s) from expectations was/were found in run {run}")
            no_deviations_in_any_run = False
    
    if no_deviations_in_any_run:
        print("No deviation from expectations was found in any of the runs")
    else:
        print("Some deviation(s) from expectations was/were found in at least one run")


if __name__ == "__main__":
    batch_csv = sys.argv[1]
    run_batch_test(batch_csv)
    
    



