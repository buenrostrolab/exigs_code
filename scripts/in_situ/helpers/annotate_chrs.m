% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Annotates match table given an ordered list of chromosomes and rules for
% unplaced/unlocalized, scaffolds

function[anno_table] = annotate_chrs(filter_table,chr_order_file)

chr_order = readtable(sprintf('%s',chr_order_file),'FileType','text');
chr_order = table2array(chr_order);

% handle unplaced/unlocalized scaffolds

filter_table.Var6(contains(filter_table.Var6,"Un_")) = {'unplaced'};
filter_table.Var6(contains(filter_table.Var6,"_random")) = {'unlocalized'};
filter_table.Var6(contains(filter_table.Var6,"_alt")) = {'unlocalized'};

% add new column of ints representing chromosome index

chr_index = zeros(length(filter_table.Var6),1);
for chr=1:length(chr_order)
   chr_index(strcmp(filter_table.Var6,chr_order(chr))) = chr;
end

filter_table.chr_idx = chr_index;
anno_table = filter_table;
anno_table = movevars(anno_table,"chr_idx",'Before',"Var6");

end