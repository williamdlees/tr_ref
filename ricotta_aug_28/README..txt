This version of the set uses revised IUIS filter criteria:

            try:
                fs_reads = float(rec['Fully_Spanning_Reads'])
                accuracy = float(rec['Percent_Accuracy'])
            except ValueError:
                continue

            # Acceptance criteria:
            # - At least 10 fully spanning reads: fs_reads >= 10
            # - Every base has at least 80% support from all spanning reads, fully matchng and otherwise: accuracy == 100
            #   ('accuracy' is is calculated as the percentage of positions in the sequence that have this level of support)
            if fs_reads >= 10 and accuracy == 100:
