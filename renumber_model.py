import argparse
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from pathlib import Path
import logging

log = logging.getLogger("log")
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)


def pdb2fastadic(pdb_f):
    """
    Write the sequence as a fasta.
    Parameters
    ----------
    pdb_f : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
    Returns
    -------
    seq_dic : dict
        dict of fasta sequences (one per chain)
    """
    res_codes = dict(
        [
            ("CYS", "C"),
            ("ASP", "D"),
            ("SER", "S"),
            ("GLN", "Q"),
            ("LYS", "K"),
            ("ILE", "I"),
            ("PRO", "P"),
            ("THR", "T"),
            ("PHE", "F"),
            ("ASN", "N"),
            ("GLY", "G"),
            ("HIS", "H"),
            ("LEU", "L"),
            ("ARG", "R"),
            ("TRP", "W"),
            ("ALA", "A"),
            ("VAL", "V"),
            ("GLU", "E"),
            ("TYR", "Y"),
            ("MET", "M"),
            ("DA", "A"),
            ("DG", "G"),
            ("DC", "C"),
            ("DT", "T"),
        ]
    )
    seq_dic = {}

    with open(pdb_f) as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                res_num = int(line[22:26])
                res_name = line[17:20].strip()
                chain = line[21]
                # if res_name in RES_TO_BE_IGNORED:
                #     continue
                try:
                    one_letter = res_codes[res_name]
                except KeyError:
                    one_letter = "X"
                if chain not in seq_dic:
                    seq_dic[chain] = {}
                seq_dic[chain][res_num] = one_letter
    return seq_dic


def align(reference, model, output_path):
    """
    Sequence align and get the numbering relationship.
    Parameters
    ----------
    reference : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
    model : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
    output_path : Path
    Returns
    -------
    align_dic : dict
        dictionary of sequence alignments (one per chain)
    """
    seqdic_ref = pdb2fastadic(reference)
    seqdic_model = pdb2fastadic(model)

    if seqdic_ref.keys() != seqdic_model.keys():
        # TODO: Implement chain-matching here
        return False

    align_dic = {}
    for ref_chain, model_chain in zip(seqdic_ref, seqdic_model):

        if ref_chain != model_chain:
            raise Exception(f"Chain mismatch: {ref_chain} != {model_chain}")

        align_dic[ref_chain] = {}

        seq_ref = Seq("".join(seqdic_ref[ref_chain].values()))
        seq_model = Seq("".join(seqdic_model[model_chain].values()))

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        alns = aligner.align(seq_ref, seq_model)
        top_aln = alns[0]

        aln_fname = Path(output_path, f"blosum62_{ref_chain}.aln")
        log.debug(f"Writing alignment to {aln_fname.name}")
        with open(aln_fname, "w") as fh:
            fh.write(str(top_aln))
        aligned_ref_segment, aligned_model_segment = top_aln.aligned

        # this should always be true
        assert len(aligned_ref_segment) == len(aligned_model_segment)

        identity = (
            str(top_aln).count("|") / float(min(len(seq_ref), len(seq_model)))
        ) * 100

        if not any(e for e in top_aln.aligned):
            # No alignment!
            log.warning(
                f"No alignment for chain {ref_chain} is it protein/dna? "
                "Matching sequentially"
            )
            if all("X" in s for s in seq_ref) and all("X" in s for s in seq_model):
                # this sequence contains only ligands, do it manually
                if len(seq_ref) != len(seq_model):
                    # we cannot handle this
                    raise f"Cannot align chain {model_chain}"
                for ref_res, model_res in zip(
                    seqdic_ref[ref_chain], seqdic_model[model_chain]
                ):

                    align_dic[ref_chain].update({model_res: ref_res})
        else:
            if identity <= 40.0:
                # Identity is very low
                log.warning(
                    f"Sequence identity of chain {ref_chain} is "
                    f"{identity:.2f}%, please check the results carefully"
                )
            else:
                log.debug(
                    f"Sequence identity between chain {ref_chain} "
                    f"of {reference} and chain {model_chain} of {model} is "
                    f"{identity:.2f}%"
                )
            for ref_segment, model_segment in zip(
                aligned_ref_segment, aligned_model_segment
            ):

                start_ref_segment, end_ref_segment = ref_segment
                start_model_segment, end_model_segment = model_segment

                reslist_ref = list(seqdic_ref[ref_chain].keys())[
                    start_ref_segment:end_ref_segment
                ]

                reslist_model = list(seqdic_model[model_chain].keys())[
                    start_model_segment:end_model_segment
                ]

                for _ref_res, _model_res in zip(reslist_ref, reslist_model):
                    align_dic[ref_chain].update({_model_res: _ref_res})

    return align_dic


def renumber_model(model_f, reference_dic):
    """
    Renumber the model according to the reference dictionary.
    """
    new_model_f = Path(model_f).name.replace(".pdb", "_renumbered.pdb")
    log.info(f"Renumbered model name: {new_model_f}")
    ignored_dic = {}
    with open(new_model_f, "w") as out_fh:
        with open(model_f, "r") as fh:
            for line in fh.readlines():
                if line.startswith("ATOM"):
                    chain = line[21]
                    resnum = int(line[22:26])
                    if chain in reference_dic:
                        if resnum in reference_dic[chain]:
                            new_resnum = reference_dic[chain][resnum]
                            line = line[:22] + str(new_resnum).rjust(4) + line[26:]
                            out_fh.write(line)
                        else:
                            if chain not in ignored_dic:
                                ignored_dic[chain] = []
                            if resnum not in ignored_dic[chain]:
                                ignored_dic[chain].append(resnum)
                else:
                    out_fh.write(line)
    for chain in ignored_dic:
        log.warning(f"Ignored residues {ignored_dic[chain]} in model's chain {chain}")
    return new_model_f


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("reference", help="")
    parser.add_argument("model", help="")
    args = parser.parse_args()
    log.setLevel("DEBUG")

    log.info("Getting sequence numbering relationship via BLOSUM62 alignment")
    numbering_dic = align(args.reference, args.model, ".")
    log.info("Renumbering model according to numbering relationship")
    renumbered_model_f = renumber_model(args.model, numbering_dic)

    log.info("Renumbering complete")
    log.info("DO NOT trust this renumbering blindly!")
    log.info("Check the .aln files for more information")
