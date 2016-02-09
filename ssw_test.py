from ssw_wrap import Aligner


polyA = Aligner("A"*200,
                match=5,
                mismatch=3,
                gap_open=4,
                gap_extend=1,
                report_secondary=False,
                report_cigar=True)


query_seq = "CTACGTAGCTAGCTAGCTATGCTAGCTGATGCTAGCTGTGTAAAAAAAAAAAAAAGAAAAAATTTAAAAAAAACGTGCTAGCTGTGCTATTAGCTAGTCGTGGCTAGTGTAC"
result = polyA.align(query_seq, min_score=20, min_len=20)
begin = result.query_begin
end = result.query_end

print(" "*begin+query_seq[begin:end+1]+"\n"+query_seq)


