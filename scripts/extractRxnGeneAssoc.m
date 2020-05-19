function [rxnGeneAssoc, rxn_groups] = extractRxnGeneAssoc(rxnAssocFile, exprDataFile, rxnScoreMethod, writeOutput)
% Extract the gene-reaction associations from the secretion model portion
%
% Input:
%
%   rxnAssocFile  File name (with path) of the .xlsx file containing the
%                 portion of the model associated with protein secretion.
%
%   writeOutput   Logical indicating whether the output (reaction-gene
%                 associations and reaction group information) should be
%                 written to .txt files. The files will be written to the
%                 same directory as in the input rxnAssocFile.
%                 (Default = True)
%
% Output:
%
%   rxnGeneAssoc  Table of reaction-gene associations, where rows are
%                 reactions, columns are genes, and a value of 1 indicates
%                 an associated gene-reaction pair and 0 otherwise.
%
%   rxn_groups    Reactions that have identical gene associations are
%                 grouped together. The rxn_groups output contains two
%                 columns, where the first column is the name of one
%                 reaction in the group, and the second column contains the
%                 names of all reactions in the group with identical gene
%                 associations.

if nargin < 2
    writeOutput = true;
end


%% Parse reaction-gene associations

% load model structure
model = importExcelModel(rxnAssocFile);

% translate gene rules from ENSG to gene abbreviations
[grRules, genes, rxnGeneMat] = translateGrRules(model.grRules, 'Names', 'ENSG');

% only keep rxns with a gene association
keep = ~cellfun(@isempty, grRules);
rxns = model.rxns(keep);
rxnGeneMat = rxnGeneMat(keep,:);

% remove reactions with identical gene associations
[~,unique_rows, rxn_groups] = unique(rxnGeneMat, 'rows', 'stable');
unique_rxns = rxns(unique_rows);
rxn_groups = arrayfun(@(i) join(rxns(i == rxn_groups), '; '), (1:numel(unique_rxns))');
rxn_groups = [unique_rxns, rxn_groups];
rxnGeneMat = rxnGeneMat(unique_rows,:);
rxns = unique_rxns;

% assemble into table
rxnGeneAssoc = array2table(rxnGeneMat, 'VariableNames', genes, 'RowNames', rxns);
rxnGeneAssoc.Properties.DimensionNames = {'Reaction', 'Gene'};

% write association table and reaction group information to file
% (writes to the same directory where the input file was located)
if ( writeOutput)
    inputPath = fileparts(rxnAssocFile);
    if isempty(inputPath)
        inputPath = '.';
    end
    writetable(rxnGeneAssoc, strcat(inputPath, '/sec_rxn_gene_assoc.txt'), 'Delimiter', '\t', 'WriteRowNames', true);
    writecell(rxn_groups, strcat(inputPath, '/sec_rxn_groups.txt'), 'Delimiter', '\t');
end


%% Compute reaction-level expression values

% load expression data file
x = readtable(exprDataFile);
x(:,1) = [];  % remove first column of indices

% extract gene names from table
% gene names are the first N column names before the "CancerStatus" column
[~,cstatus_indx] = ismember('CancerStatus', x.Properties.VariableNames);
expr_genes = x.Properties.VariableNames(1:(cstatus_indx - 1))';
expr_levels = table2array(x(:,1:(cstatus_indx - 1)));

% calculate reaction scores from gene expression values
nrow = size(x,1);
rxnScores = zeros(nrow, numel(rxns));  % initialize
switch rxnScoreMethod
    case 'mean'
        for i = 1:nrow
            rxnScores(i,:) = arrayfun(@(a) mean(expr_levels(i, rxnGeneMat(a,:)==1)), 1:numel(rxns));
        end
end














