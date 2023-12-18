import os
import re
from tqdm import tqdm
from toolbiox.lib.common.util import printer_list
from toolbiox.api.common.mapping.blast import outfmt5_read_big, keep_outfmt6_info, outfmt5_complete
from toolbiox.lib.common.fileIO import tsv_file_dict_parse
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Gff2Gtf(object):
    u"""
    convert gff3 to gtf
    __author__ = "Zhang Yiming"
    """

    def __init__(self, input, output):
        self.input = input
        self.output = output

        self.genes = {}
        self.transcripts = {}
        self.convert()
        pass

    @staticmethod
    def split_gff_detail(line):
        u"""
        将gff中的信息列，拆分成字典
        """
        res = {}

        for i in line.split(";"):
            i = i.split("=")

            res[i[0]] = i[1] if ":" not in i[1] else i[1].split(":")[1]
        return res

    @staticmethod
    def concat_dict_to_string(data):
        u"""
        将字典中的键值对构成字符串
        """
        res = []

        for k, v in data.items():
            if "id" in k or "name" in k:
                res.append("%s \"%s\"" % (k, v))
        return "; ".join(res)

    def convert(self):
        u"""
        start to converting
        """
        with open(self.output, "w+") as w:
            w.write("#gtf-version")
            with open(self.input) as r:
                for line in tqdm(r):
                    # skip empty line
                    if not line.strip():
                        continue

                    if line.startswith("#"):
                        w.write(line)
                        continue

                    lines = line.rstrip().split("\t")
                    info = self.split_gff_detail(lines[-1])

                    # first class. eg: gene
                    if "ID" in info.keys() and \
                        "Parent" not in info.keys() and \
                            "gene" in lines[2]:

                        if "gene_id" not in info.keys():
                            info["gene_id"] = info["ID"]

                        if "gene_name" not in info.keys():
                            info["gene_name"] = info["Name"]

                        self.genes[info["gene_id"]] = info["gene_name"]

                    elif "Parent" in info.keys():

                        # second class. eg: transcripts
                        if info["Parent"] in self.genes.keys():
                            info["gene_id"] = info["Parent"]
                            info["gene_name"] = self.genes[
                                info["Parent"]
                            ]

                            if "transcript_id" not in info.keys():
                                info["transcript_id"] = info["ID"]

                            if "transcript_name" not in info.keys():
                                info["transcript_name"] = info["Name"] if "Name" in info.keys(
                                ) else "NA"

                            self.transcripts[info["transcript_id"]] = [
                                info["transcript_name"],
                                info["gene_id"]
                            ]

                            info["transcript_type"] = lines[2]

                            if lines[2] != "CDS":
                                lines[2] = "transcript"

                        # third class. eg: exons
                        else:
                            if "transcript_id" not in info.keys():
                                info["transcript_id"] = info["Parent"]

                            transcript = self.transcripts[
                                info["Parent"]
                            ]

                            if "transcript_name" not in info.keys():
                                info["transcript_name"] = transcript[0]

                            if "gene_id" not in info.keys():
                                info["gene_id"] = transcript[1]

                            if "gene_name" not in info.keys():
                                info["gene_name"] = self.genes[
                                    transcript[1]
                                ]

                            ele_id = "%s_id" % lines[2].lower()
                            ele_name = "%s_name" % lines[2].lower()

                            if ele_id not in info.keys() and \
                                    "ID" in info.keys():
                                info[ele_id] = info.pop("ID")

                            if ele_name not in info.keys() and \
                                    "Name" in info.keys():
                                info[ele_name] = info.pop("Name")

                    # others just convert and output
                    else:
                        pass

                    lines[-1] = self.concat_dict_to_string(info)

                    if lines[2] in ("gene", "transcript", "exon", "CDS"):
                        w.write("\t".join(lines) + "\n")


class Gtf2Gff(object):
    """
    __author__ = "Zhang Yiming"
    """

    def __init__(self, input, output, generate_genes=False):
        u"""
        init this class
        """
        self.input = input
        self.output = output
        self.generate_genes = generate_genes

        self.genes = {}
        self.transcripts = {}
        self.convert()
        pass

    @staticmethod
    def __split_gtf_details__(line):
        u"""
        format gtf detailed messages into dict
        :param line: columns like gene_id "gene"; transcript_id "transcript";
        :return: dict {"gene_id": "gene"}
        """
        data = {}
        for message in line.strip().split(";"):
            if not message:
                continue
            key, value = message.strip().split(" ")
            data[key.lower()] = re.sub(r"[\";]", "", value)
        return data

    @staticmethod
    def __get_value_from_data__(data, target, pop=True):
        u"""
        从gtf的详细列中构建出的字典，从其中提取出所需要的数据
        但是由于有多个可能性，比如：gene_id, geneID, ID等等，不同的标准下的gtf文件，太烦人了
        因此，在此通过正则来处理这个问题
        :param data: 从gtf文件中，提取出的字典信息
        :param target: 所要提取数据的目标，为gene_id等标签
        :return: string
        """
        ids = [target, target.replace("_", ""), target.split("_")[-1]]
        for i in ids:
            if i in data.keys():
                return data.pop(i) if pop else data[i]
        return "NA"

    def __format_gff_details__(self, data, label):
        u"""
        将获取到的gtf的信息，format成gff3样式
        :param data: 由self.__split_gtf_details__构造的字典
        :param label: gtf文件，第二列表明的元件类型
        :return: string
        """
        result = ""
        if label == "gene":
            result += "ID=%s;Name=%s" % (
                self.__get_value_from_data__(data=data, target="gene_id"),
                self.__get_value_from_data__(data=data, target="gene_name"),
            )
        elif label == "transcript":
            result += "ID=%s;Name=%s;Parent=%s" % (
                self.__get_value_from_data__(
                    data=data, target="transcript_id"),
                self.__get_value_from_data__(
                    data=data, target="transcript_name"),
                self.__get_value_from_data__(data=data, target="gene_id"),
            )

            if "gene_name" in data.keys():
                data.pop("gene_name")
        else:
            ids = ["exon_id", "protein_id", "%s_id" % label.lower(), "ID"]

            for i in ids:
                if i in data.keys():
                    ids = data.pop(i)
                    break

            ids = ids if isinstance(ids, str) else "NA"

            parent = self.__get_value_from_data__(
                data=data, target="transcript_id")

            if "exon_number" in data.keys() and ids == "NA":
                ids = "%s.%s" % (parent, data["exon_number"])

            result += "ID=%s;Parent=%s" % (ids, parent)

            if "gene_name" in data.keys():
                data.pop("gene_name")

        for key, value in data.items():
            result += ";%s=%s" % (key, value)

        return result

    def convert(self):
        u"""
        进行转化
        :return:
        """
        genes = set()
        with open(self.output, "w+") as w:
            w.write("#gff-version 3\n")
            with open(self.input) as r:
                for line in r:
                    if line.startswith("#"):
                        continue

                    # skip empty line
                    if not line.strip():
                        continue

                    lines = line.split("\t")
                    data = self.__split_gtf_details__(lines[8])

                    if self.generate_genes:
                        if lines[2] == "gene":
                            genes.add(self.__get_value_from_data__(
                                data, "gene_id", False))

                        if lines[2] == "transcript":
                            parent = self.__get_value_from_data__(
                                data, "gene_id", False)
                            parent_name = self.__get_value_from_data__(
                                data, "gene_name", False)

                            if parent_name == "NA":
                                parent_name = parent
                            if parent not in genes:
                                new_line = lines[:8] + ["ID=%s;Name=%s" %
                                                        (parent, parent_name)]
                                new_line[2] = "gene"
                                genes.add(parent)
                                w.write("\t".join(new_line) + "\n")

                    lines[8] = self.__format_gff_details__(data, lines[2])

                    w.write("\t".join(lines) + "\n")


def fancy_name_parse(input_string):
    contig_name, c_start, c_end = re.search(
        r'^(\S+):(\d+)\.\.(\d+)$', input_string).groups()
    return contig_name, int(c_start), int(c_end)


def genblasta2BED_main(args):
    """
    class abc(object):
        pass

    args = abc()

    args.input_file = '/lustre/home/xuyuxing/Work/Csp/ITS/Cau.rRNA'
    args.output_file = '/lustre/home/xuyuxing/Work/Csp/ITS/Cau.rRNA.bed'
    args.ID_prefix = 'Cau_ITS_'
    """

    gb_file = tsv_file_dict_parse(args.input_file, seq="|",
                                  fieldnames=['query_name', 'subject_name', 'strand', 'gene_cover', 'score',
                                              'rank'])

    with open(args.output_file, 'w') as f:
        num = 0
        for i in gb_file:
            if gb_file[i]['rank'] is not None and re.match(r'rank:\d+', gb_file[i]['rank']):
                num = num + 1

                score = float(
                    re.search(r'score:(.*)', gb_file[i]['score']).group(1))

                contig_name, c_start, c_end = fancy_name_parse(
                    gb_file[i]['subject_name'])

                f.write("%s\t%d\t%d\t%s\t%f\t%s\n" % (
                    contig_name, c_start, c_end, args.ID_prefix + str(num), score, gb_file[i]['strand']))


class BlastFormat(object):
    """
    blast格式转换
    """

    outfmt6_fieldnames = ["query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_openings",
                          "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"]

    def __init__(self, input_blast_file=None, input_fasta_file=None, output_file=None):
        self.input_blast_file = input_blast_file
        self.input_fasta_file = input_fasta_file
        self.output_file = output_file

    def outfmt5To6(self):
        """
        blast outfmt 5 to 6
        :return:
        """
        output_dict = outfmt5_read_big(self.input_blast_file, False)
        with open(self.input_blast_file, 'w') as f:
            for query in output_dict:
                for hsp in keep_outfmt6_info(query):
                    f.write(printer_list(hsp) + "\n")

    def outfmt5complete(self):
        """
        判断blast outfmt 5文件是否完整
        :return:
        """

        if outfmt5_complete(self.input_blast_file):
            print("%s is complete" % self.input_blast_file)
        else:
            print("%s is not complete" % self.input_blast_file)

    def outfmt6ToFasta(self):
        """
        blast outfmt 6 to fasta
        :return:
        """

        blast_file = tsv_file_dict_parse(
            self.input_blast_file, fieldnames=self.outfmt6_fieldnames)

        ref_dict = Fasta(self.input_fasta_file)

        with open(self.output_file, 'w') as f:
            for ID in blast_file:
                s_name = blast_file[ID]['subject_id']
                s_start = int(blast_file[ID]['s_start'])
                s_end = int(blast_file[ID]['s_end'])

                if s_end > s_start:
                    neg_strand = False
                    strand = "+"
                else:
                    neg_strand = True
                    strand = "-"
                    tmp = s_start
                    s_start = s_end
                    s_end = tmp

                a = ref_dict.get_seq(s_name, s_start, s_end, rc=neg_strand)
                fancy_name = "%s:%d-%d:%s" % (s_name, s_start, s_end, strand)

                contig_record = SeqRecord(
                    Seq(a.seq), id=ID, description=fancy_name)

                f.write(contig_record.format("fasta"))
