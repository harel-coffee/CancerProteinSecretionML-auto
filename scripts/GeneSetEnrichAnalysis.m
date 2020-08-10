function gsea_results = GeneSetEnrichAnalysis(gene_list,gene_sets,all_genes,gs_bounds,depletion)
% Performs a Gene Set Enrichment Analysis
%
%--------------------------------- INPUTS ---------------------------------
%
% gene_list     list of genes that will be evaluated for enrichment of 
%               certain gene sets.
%
% gene_sets     cell array defining the gene sets, where the first column
%               contains the name of the gene sets, and the second column
%               contains the names of the genes belonging to each set. 
%
%                   e.g.: GENE_SETS = { 'gene set 1', 'GENE A'
%                                       'gene set 1', 'GENE B'
%                                       'gene set 1', 'GENE C'
%                                       'gene set 2', 'GENE A'
%                                       'gene set 2', 'GENE D'
%                                       'gene set 3', 'GENE E'
%                                       'gene set 3', 'GENE F' }
%
%           *** NOTE: rows of GENE_SETS that do not contain a gene from
%               GENE_LIST or ALL_GENES will be removed before the analysis.
%               This does not impact the calculations, but will affect the
%               "set size" value in the GSEA_RESULTS structure.
%
% all_genes     [OPTIONAL] a list of all possible genes from which 
%               GENE_LIST and GENE_SETS are drawn. If excluded, ALL_GENES
%               will be set as the list of all genes contained in GENE_LIST
%               and GENE_SETS.
%
% gs_bounds     (Default = no bounds)
%
% depletion     (Default = False) If TRUE, instead of testing the list of
%               genes for ENRICHMENT of certain sets, it will instead test
%               the list for DEPLETION of certain sets.
%
%--------------------------------- OUTPUTS --------------------------------
%
% gsea_results  a cell array containing the results of the enrichment
%               analysis, with gene sets sorted by ascending p-value.
%


%% Handle input arguments
if nargin < 5
    depletion = false;
    if nargin < 4 || isempty(gs_bounds)
        gs_bounds = [0 Inf];
        if nargin < 3 || isempty(all_genes)
            % ALL_GENES will be the combination of all genes in GENE_LIST and GENE_SETS
            all_genes = unique([gene_list;gene_sets(:,2)]);
        end
    end
end

% set min and max gene set (GS) sizes
[minGSsize,maxGSsize] = deal(gs_bounds(1),gs_bounds(2));


%% Validate inputs

% check for duplicated gene names
if length(unique(gene_list)) < length(gene_list)
    error('Duplicate entries exist in GENE_LIST.');
end

% check for GENE_LIST entries not included in ALL_GENES
if any(~ismember(gene_list,all_genes))
    error('GENE_LIST contains genes that are not present in ALL_GENES.');
end


%% Pre-process gene sets (GS)

% remove rows of GENE_SETS that do not contain genes in ALL_GENES
fprintf('Checking for empty gene sets... ');
a = length(unique(gene_sets(:,1)));  % check number of sets before removal
ind = ~ismember(gene_sets(:,2),all_genes);
gene_sets(ind,:) = [];
b = length(unique(gene_sets(:,1)));  % check number of sets after removal
fprintf('Removed %u empty sets.\n',a-b);

% remove repeated rows in GENE_SETS
fprintf('Checking for duplicated rows in GENE_SETS... ');
a = size(gene_sets,1);  % check number of rows before removal
[~,gs_ind] = ismember(gene_sets,gene_sets);
[~,uniq_ind] = unique(gs_ind,'rows');
gene_sets = gene_sets(uniq_ind,:);
b = size(gene_sets,1);  % check number of rows after removal
fprintf('Removed %u duplicated rows.\n',a-b);

% determine gene set sizes and remove those not satisfying constraints
fprintf('Checking gene set sizes... ');
[GSnames,GSsizes] = cellfreq(gene_sets(:,1));
ind = (GSsizes < minGSsize) | (GSsizes > maxGSsize);
gene_sets(ismember(gene_sets(:,1),GSnames(ind)),:) = [];
GSnames(ind) = []; GSsizes(ind) = [];
fprintf('Removed %u gene sets not satisfying size limits.\n',sum(ind));

% convert GENES and GENE_SETS to numeric arrays to speed up calculations
% gene_sets_proc = gene_sets;  % save processed GENE_SETS
% GSnames = unique(gene_sets(:,1),'stable');  % save list of gene set names

[~,set_nums] = ismember(gene_sets(:,1),gene_sets(:,1));  % index gene sets
[~,g_nums] = ismember(gene_sets(:,2),all_genes);  % index genes
gene_sets_num = [set_nums,g_nums];  % numeric GENE_SETS
[~,GSnums] = ismember(GSnames,gene_sets(:,1));  % unique list of gene set numbers

[~,gene_list_num] = ismember(gene_list,all_genes);  % numeric GENE_LIST
[~,all_genes_num] = ismember(all_genes,all_genes);  % numeric ALL_GENES

% [GSnums,GSsizes] = cellfreq(gene_sets_num(:,1));  % recalc gene set sizes


%% Calculate significance of gene set enrichments

fprintf('Calculating enrichment... ');
[penr,pdep] = arrayfun( @(g) EnrichmentTest(all_genes_num,gene_list_num,gene_sets_num(ismember(gene_sets_num(:,1),g),2)),GSnums);
n_hits = arrayfun( @(g) length(intersect(gene_list_num,gene_sets_num(ismember(gene_sets_num(:,1),g),2))),GSnums);

if ( depletion )
    p = pdep;
else
    p = penr;
end
fprintf('Done.\n');


%% Process and organize results

% sort results in order of ascending p-value
[p,sort_ind] = sort(p);
GSsizes = GSsizes(sort_ind);
n_hits = n_hits(sort_ind);
GSnames = GSnames(sort_ind);

% calculate adjusted p-values
p_holm = adjust_pvalues(p,'holm');
p_bh = adjust_pvalues(p,'benjamini');

% collect results into cell array
gsea_results = [{'set name','set size','no. hits','p (raw)','padj (Holm)','padj (BH)'}; ...
    GSnames, num2cell([GSsizes, n_hits, p, p_holm, p_bh])];


end


%% Additional functions

function [item,freq] = cellfreq(c)
%Calculate the frequency of each item in a cell array.
item = unique(c,'stable');  % determine unique items in C
[~,c] = ismember(c,item);  % convert C to numeric values to increase speed
edges = 0.5:1:(length(item)+0.5);  % define edges between each value of C
freq = histcounts(c,edges)';  % count occurrence of each value of C
end

