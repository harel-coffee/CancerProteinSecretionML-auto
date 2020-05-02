% Load and analyze data/results from Cancer Omics Data Exploration project


%% Load ranking/score data for separating cancers by type (normal v. tumor)


%********************* OPTIONS *********************
group_by = 'model';  % 'cancer' or 'model'
filter_by = 'rocauc';  % 'rocauc' or 'accuracy'
filter_models = {'ExtraTreesClassifier';'RandomForestClassifier';'AdaBoostClassifier';'XGBClassifier';'SVC';'LogisticRegression'};  % list of models to include in the average model score used for filtering (leave empty to include all)
filter_thresh = 0.75; %0.75;  % average model score necessary to include cancer
scale_by_rocauc = false;  % if TRUE, then gene-level scores will be scaled by the ROC AUC
                         % according to the following formula:
                         % scaled_score = score * (2*ROCAUC - 1)

select_dataset = 'mutTP53';
% WARNING: Case-sensitive!
%
%   'CancerStatus': cancer status (normal vs. tumor)
%
%   'mutGENE': compare tumor samples with/without GENE mutation
%              (e.g., 'mutTP53')
%
%   'allmuts': load all available gene mutation comparisons
%
%   'stageXvY': compare tumor stage X with stage Y
%               Options =  'stage1v2', 'stage1v3', 'stage1v4',
%                          'stage2v3', 'stage2v4', 'stage3v4',
%                          'stage1v234', 'stage12v34', 'stage123v4'
%               Note: 'stage1v234' indicates that stage 1 was compared to a
%                     group of samples that included stages 2, 3, and 4.
%
%   'allstages': load all available tumor stage comparisons

exclude_cancers = {};
switch lower(select_dataset)
    case 'allmuts'
        file_phrase = 'mut*';
    case 'stage1v2'
        file_phrase = 'TumorStage_stagei_stageii';
    case 'stage1v3'
        file_phrase = 'TumorStage_stagei_stageiii';
    case 'stage1v4'
        file_phrase = 'TumorStage_stagei_stageiv';
    case 'stage2v3'
        file_phrase = 'TumorStage_stageii_stageiii';
    case 'stage2v4'
        file_phrase = 'TumorStage_stageii_stageiv';
    case 'stage3v4'
        file_phrase = 'TumorStage_stageiii_stageiv';
    case 'stage1v234'
        file_phrase = 'TumorStage_stagei_stageii-iv';
    case 'stage12v34'
        file_phrase = 'TumorStage_stagei-ii_stageiii-iv';
    case 'stage123v4'
        file_phrase = 'TumorStage_stagei-iii_stageiv';
    case 'allstages'
        file_phrase = 'TumorStage_*';
    otherwise
        file_phrase = select_dataset;
end


% get file information from selected directory
fileinfo = dir(['/Users/jonrob/Documents/PostDoc/CancerOmicsDataExploration/results/*/*_',file_phrase,'_*']);

% initialize variables
[xdata,xtext,de,genes,cancers,alldata,modelscores] = deal({});

% extract data from each file
modeltypes = {};
fprintf('Extracting data from files... ');
for i = 1:length(fileinfo)
    
    % get cancer type from file name
    if strcmp(select_dataset, 'allmuts')
        % include mutated gene name in cancer name
        tmp = strsplit(fileinfo(i).name,'_');
        cname = {[tmp{1},'-',tmp{2}]};
        
    elseif strcmp(select_dataset, 'allstages')
        % shorten stage labels to numbers and append to cancer name
        stage_rename = {'stagei',       '1';
                        'stageii',      '2';
                        'stageiii',     '3';
                        'stageiv',      '4';
                        'stagei-ii',    '12';
                        'stagei-iii',   '123';
                        'stageii-iv',   '234';
                        'stageiii-iv',  '34'};
        tmp = strsplit(fileinfo(i).name,'_');
        [hasMatch,ind] = ismember(tmp, stage_rename(:,1));
        tmp(hasMatch) = stage_rename(ind(hasMatch),2);
        cname = {[tmp{1}, '-stage', [tmp{3},'v', tmp{4}]]};
    else
        % include only cancer type
        tmp = strsplit(fileinfo(i).name,'_');
        cname = tmp(1);
    end
    
    % determine or assign cancer index
    [~,cancer_ind] = ismember(cname,cancers);
    if cancer_ind == 0
        cancers = [cancers; cname];
        cancer_ind = length(cancers);
    end
    
    if strfind(fileinfo(i).name,'GenesRanking')
        
        % import gene ranking data from file
        x = importdata([fileinfo(i).folder,'/',fileinfo(i).name]);
        xdata{cancer_ind,1} = x.data;
        xtext{cancer_ind,1} = x.textdata;
        
        % obtain list of model types corresponding to gene scores
        modeltypes = x.textdata(1,2:end)';
        
        % assemble list of all genes among all data files
        genes = union(genes, x.textdata(2:end,1));
        
    elseif strfind(fileinfo(i).name,'Accuracy')
        
        % import model accuracy data from file
        x = importdata([fileinfo(i).folder,'/',fileinfo(i).name]);
        modelscores.accuracy(:,cancer_ind) = x.data(:,1);
        modelscores.accuracyCIlower(:,cancer_ind) = x.data(:,2);
        modelscores.accuracyCIupper(:,cancer_ind) = x.data(:,3);
        
    elseif strfind(fileinfo(i).name,'ROCAUC')
        
        % import model ROC AUC data from file
        x = importdata([fileinfo(i).folder,'/',fileinfo(i).name]);
        modelscores.rocauc(:,cancer_ind) = x.data(:,1);
        modelscores.rocaucCIlower(:,cancer_ind) = x.data(:,2);
        modelscores.rocaucCIupper(:,cancer_ind) = x.data(:,3);
        
        % get list of model types corresponding to model scores
        if ~isfield(modelscores,'names')
            % only have to assign once, since they should be identical in every file
            modelscores.names = x.textdata(2:end,1);
        end
        
    elseif strfind(fileinfo(i).name,'DEresults')
        
        % import differential expression (DE) analysis results from file
        x = readtable([fileinfo(i).folder,'/',fileinfo(i).name]);
        de.data{cancer_ind,1} = [x.logFC, x.FDR];
        de.datatypes = {'DElogfc';'DEfdr'};
        de.genes{cancer_ind,1} = x.pspGenes;
    
        % add to list of all genes among all data files
        genes = union(genes, de.genes{cancer_ind});
        
    end
    
end
fprintf('Done.\n');
alldata.genes = genes;

% identify cancers with poor model scores
if ~isempty(modelscores)
    
    if isempty(filter_models)
        model_ind = true(length(modelscores.names),1);
    else
        model_ind = ismember(modelscores.names,filter_models);
    end
    
    if strcmpi(filter_by,'accuracy')
        meanscores = mean(modelscores.accuracy(model_ind,:))';
    elseif strcmpi(filter_by,'rocauc')
        meanscores = mean(modelscores.rocauc(model_ind,:))';
    end
    
    % add poor-scoring cancers to 'exclude' list
    exclude_cancers = unique([exclude_cancers; cancers(meanscores < filter_thresh)]);
    rem_ind = ismember(cancers,exclude_cancers);
    
    % remove score data associated with poor-scoring cancers
    scorefields = fields(modelscores);
    scorefields(ismember(scorefields,'names')) = [];  % ignore "names" field
    for i = 1:length(scorefields)
        modelscores.(scorefields{i})(:,rem_ind) = [];
    end
    
end

% remove remaining data for excluded cancers
rem_ind = ismember(cancers,exclude_cancers);
cancers(rem_ind) = [];
xdata(rem_ind) = [];
xtext(rem_ind) = [];
if ~isempty(de)
    de.data(rem_ind) = [];
    de.genes(rem_ind) = [];
else
    de.data = [];
    de.datatypes = [];
    de.genes = {};
end

% scale gene scores by ROC AUC if requested
if scale_by_rocauc
    for i = 1:numel(cancers)
        
        % calculate scale factors for each model
        scale_factors = 2*(modelscores.rocauc(:,i)') - 1;
        
        % scale gene scores
        [~,model_ind] = ismember(modelscores.names, modeltypes);
        gene_scores = xdata{i}(:,model_ind) .* scale_factors;
        xdata{i}(:,model_ind) = gene_scores;
        
        % re-calculate average model scores
        xdata{i}(:,end) = mean(gene_scores,2);
        
    end
end


% assemble data structure
fprintf('Organizing data... ');
if strcmp(group_by,'cancer')
    
    % add list of model types to data structure
    alldata.modeltypes = [modeltypes; de.datatypes];
    
    % add data to each field
    for i = 1:numel(cancers)
        
        % initialize structure fields with zeros
        alldata.(cancers{i}) = zeros(length(genes),numel(alldata.modeltypes));
        alldata.(cancers{i})(:,ismember(alldata.modeltypes,'DEfdr')) = 1;  % initialize p-values with 1 instead of 0
        
        % add results from machine learning models
        [~,g_ind] = ismember(xtext{i}(2:end,1),alldata.genes);
        [~,m_ind] = ismember(xtext{i}(1,2:end),alldata.modeltypes);
        alldata.(cancers{i})(g_ind,m_ind) = xdata{i};
        
        % add results from DE analysis
        [~,g_ind] = ismember(de.genes{i},alldata.genes);
        [~,m_ind] = ismember(de.datatypes,alldata.modeltypes);
        alldata.(cancers{i})(g_ind,m_ind) = de.data{i};
        
    end
    
    % add normalized, log-transformed FDR values as a new score type
    [~,fdr_ind] = ismember('DEfdr',alldata.modeltypes);
    if fdr_ind ~= 0
        alldata.modeltypes(end+1) = {'log10DEfdr'};
        for i = 1:numel(cancers)
            x = -log10(alldata.(cancers{i})(:,fdr_ind));  % log-transform
            alldata.(cancers{i})(:,end+1) = x./max(x,[],1);  % normalize by max value in each cancer type
        end
    end
    
    % calculate average gene scores among some of the higher-scoring model types
    alldata.modeltypes(end+1) = {'goodavg'};
    goodmodels = {'ExtraTreesClassifier';'RandomForestClassifier';'AdaBoostClassifier';
                  'XGBClassifier';'SVC';'LogisticRegression'};
    goodmodel_inds = ismember(modeltypes,goodmodels);
    for i = 1:numel(cancers)
        alldata.(cancers{i})(:,end+1) = mean(alldata.(cancers{i})(:,goodmodel_inds),2);
    end
    
elseif strcmp(group_by,'model')
    
    % add list of cancers to data structure
    alldata.cancers = cancers;
    
    % add DE methods to list of modeltypes
    modeltypes = [modeltypes; de.datatypes];
    
    % initialize structure fields with zeros
    for i = 1:length(modeltypes)
        alldata.(modeltypes{i}) = zeros(length(genes),length(cancers));
    end
    
    % add data to each field
    for i = 1:length(alldata.cancers)
        
        [~,g_ind_ML] = ismember(xtext{i}(2:end,1),alldata.genes);
        mtypes = xtext{i}(1,2:end)';
        for j = 1:length(mtypes)
            alldata.(mtypes{j})(g_ind_ML,i) = xdata{i}(:,j);
        end
        
        if (~isempty(de.data)) && (numel(de.genes) >= i) && (~isempty(de.genes{i}))
            [~,g_ind_DE] = ismember(de.genes{i},alldata.genes);
            for j = 1:length(de.datatypes)
                alldata.(de.datatypes{j})(g_ind_DE,i) = de.data{i}(:,j);
            end
        end
        
    end
    
    % add normalized log-transformed FDR values as a new field
    if isfield(alldata,'DEfdr')
        x = -log10(alldata.DEfdr);  % log-transform
        alldata.log10DEfdr = x./max(x,[],1);  % normalize by max value in each cancer type
    end
    
    % calculate average gene scores among some of the higher-scoring model types
    goodmodels = {'ExtraTreesClassifier';'RandomForestClassifier';'AdaBoostClassifier';
                  'XGBClassifier';'SVC';'LogisticRegression'};
    goodavg = [];
    for i = 1:length(goodmodels)
        goodavg(:,:,i) = alldata.(goodmodels{i});
    end
    alldata.goodavg = mean(goodavg,3);
    
end
fprintf('Done.\n');

% add model scores to alldata structure
alldata.modelscores = modelscores;


% clear intermediate variables
clearvars -except alldata


%% Evaluate correlation of gene scores between the different methods

% specify fieldnames
f = {'ExtraTreesClassifier';'RandomForestClassifier';'AdaBoostClassifier';
    'XGBClassifier';'LinearDiscriminantAnalysis';'SVC';'LogisticRegression';
    'log10DEfdr'};

% iterate through cancer types
C = zeros(numel(f), numel(f), numel(alldata.cancers));  % initialize
P = C;
for i = 1:numel(alldata.cancers)
    
    % organize data into matrix
    sel_data = zeros(numel(alldata.genes), numel(f));
    for j = 1:numel(f)
        sel_data(:,j) = alldata.(f{j})(:,i);
        if strcmpi(f{j},'DElogfc')
            sel_data(:,j) = abs(sel_data(:,j));  % take abs val of fold-changes
        end
    end
    
    % calculate correlation
    [c,p] = corr(sel_data, 'type', 'Spearman');
    C(:,:,i) = c;
    P(:,:,i) = p;
    
end

% visualize with heatmap
cmap = custom_cmap('whitemagma');

% single cancer type
cancerType = 'LGG';
genHeatMap(C(:,:,ismember(alldata.cancers,cancerType)),f,f,'both','euclid',cmap,[0 1],'k');

% mean over all cancer types
Cavg = mean(c,3);
Cmed = median(c,3);
genHeatMap(Cmed,f,f,'both','euclid',cmap,[0 1],'k');


%------- Compare cancers and models using dimensionality reduction --------

% collect data
label_cancer = {};
label_model = {};
mergedData = [];
for i = 1:numel(alldata.cancers)
    for j = 1:numel(f)
        label_cancer = [label_cancer; alldata.cancers{i}];
        label_model = [label_model; f{j}];
        if strcmpi(f{j},'DElogfc')
            x = abs(alldata.(f{j})(:,i));
            x = x./max(x);  % normalize values to 0-1 scale
            mergedData = [mergedData, x];
        else
            mergedData = [mergedData, alldata.(f{j})(:,i)];
        end
    end
end

% run PCA
[pca_coeff,pca_score,~,~,pca_pctvar] = pca(mergedData','NumComponents',3);

% run tSNE
tsne_coords = tsne(mergedData', 'NumDimensions', 2, 'Perplexity', 10);

% plot results
labels = label_cancer;
grouped_scatter(pca_score, labels, 3, 'lines', 50);
set(gca, 'FontSize', 12);
legend(unique(labels), 'FontSize', 12);


%% Compare model ROC AUC among different cancer or model types

% generate boxplots (group by cancer type)
figure;
[~,sort_ind] = sort(-median(alldata.modelscores.rocauc,1));
if size(alldata.modelscores.rocauc, 2) > 20
    boxplot(alldata.modelscores.rocauc(:,sort_ind), 'PlotStyle', 'compact', 'Colors', 'k');
else
    boxplot(alldata.modelscores.rocauc(:,sort_ind), 'Symbol', 'ko', 'Jitter', 0.25);
end
set(gca,'XTick',1:numel(alldata.cancers),'XTickLabel',alldata.cancers(sort_ind),'XTickLabelRotation',90);

% generate boxplots (group by model type)
figure;
[~,sort_ind] = sort(-median(alldata.modelscores.rocauc,2));
if size(alldata.modelscores.rocauc, 1) > 20
    boxplot(alldata.modelscores.rocauc(sort_ind,:)', 'PlotStyle', 'compact', 'Colors', 'k');
else
    boxplot(alldata.modelscores.rocauc(sort_ind,:)', 'Symbol', 'ko', 'Jitter', 0.25);
end
set(gca,'XTickLabel',alldata.modelscores.names(sort_ind),'XTickLabelRotation',90);

% generate heatmap
figure;
cmap = flipud(custom_cmap('whitemagma'));
genHeatMap(alldata.modelscores.rocauc', alldata.modelscores.names, alldata.cancers, 'both', 'euclid', cmap, [], 'k');


%% Visualize top-scoring genes across cancer types

% specify parameters
min_gene_score = 0.5;  % min gene score to include gene and/or cancer
score_method = 'min';  % 'min', 'mean', or 'median'
model_type = 'Average';  % name of model (e.g., 'XGBClassifier'), or 'Average'
cmap = flipud(custom_cmap('whitemagma'));  % colormap

% filter genes and cancer types
cind = true(size(alldata.cancers));
if strcmpi(score_method,'min')
    cind = any(alldata.(model_type) > min_gene_score, 1);  % optional cancer type filtration
    gind = any(alldata.(model_type) > min_gene_score, 2);
elseif strcmpi(score_method,'mean')
    gind = mean(alldata.(model_type), 2) > min_gene_score;
elseif strcmpi(score_method,'median')
    gind = median(alldata.(model_type), 2) > min_gene_score;
end


% generate heatmap
figure
genHeatMap(alldata.(model_type)(gind,cind), alldata.cancers(cind), alldata.genes(gind), 'both', 'euclid', cmap, [], 'k');


% generate boxplot
[~,sort_gene_ind] = sort(median(alldata.(model_type)(gind,cind),2));
plotdata = alldata.(model_type)(gind,cind)';
plotdata = plotdata(:, sort_gene_ind);
genenames = alldata.genes(gind);
genenames = genenames(sort_gene_ind);
figure
boxplot(plotdata, 'PlotStyle', 'compact', 'Orientation', 'horizontal', 'Colors', 'k');
set(gca,'YTickLabel', genenames, 'YTick', 1:numel(genenames));%,'XTickLabelRotation',90);


%% Visualize top-scoring genes across models within ONE CANCER type

% specify parameters
min_gene_score = 0.5;  % min gene score to include gene
score_method = 'min';  % 'min', 'mean', or 'median'
cancer_type = 'COAD';
f = {'ExtraTreesClassifier';'RandomForestClassifier';'AdaBoostClassifier';
     'XGBClassifier';'LinearDiscriminantAnalysis';'SVC';'LogisticRegression';
     'Average';'log10DEfdr'};

% collect data for selected cancer type
[~,cancer_ind] = ismember(cancer_type, alldata.cancers);
sel_data = [];
for i = 1:numel(f)
    sel_data(:,i) = alldata.(f{i})(:,cancer_ind);
end

% filter genes
if strcmpi(score_method,'min')
    gind = any(sel_data > min_gene_score, 2);
elseif strcmpi(score_method,'mean')
    gind = mean(sel_data, 2) > min_gene_score;
elseif strcmpi(score_method,'median')
    gind = median(sel_data, 2) > min_gene_score;
end

% generate heatmap
figure
cmap = flipud(custom_cmap('whitemagma'));
genHeatMap(sel_data(gind,:), f, alldata.genes(gind), 'both', 'euclid', cmap, [], 'k');


% generate boxplot
[~,sort_gene_ind] = sort(median(sel_data(gind,:),2));
plotdata = sel_data(gind,:)';
plotdata = plotdata(:, sort_gene_ind);
genenames = alldata.genes(gind);
genenames = genenames(sort_gene_ind);
figure
boxplot(plotdata, 'PlotStyle', 'compact', 'Orientation', 'horizontal', 'Colors', 'k');
set(gca,'YTickLabel', genenames, 'YTick', 1:numel(genenames));


%% Visualize scores for ONE GENE across models and cancer types

% specify gene name and model types to examine
gene = 'HSPA4L';  % specify gene name
f = {'ExtraTreesClassifier';'RandomForestClassifier';'AdaBoostClassifier';
    'XGBClassifier';'LinearDiscriminantAnalysis';'SVC';'LogisticRegression';
    'log10DEfdr'};

% collect gene scores across cancer and model types
[~,g_ind] = ismember(gene, alldata.genes);
gene_scores = [];
for i = 1:numel(f)
    gene_scores(:,i) = alldata.(f{i})(g_ind,:);
end

% generate heatmap
cmap = flipud(custom_cmap('whitemagma'));
genHeatMap(gene_scores,f,alldata.cancers,'both','euclid',cmap,[],'k');


%% Visualize top-scoring genes for ONE MODEL across cancer types

model = 'log10DEfdr';  % specify model name
min_gene_score = 0.8;
score_method = 'min';  % 'min', 'mean', or 'median'

sel_data = alldata.(model);

% filter genes
cind = true(size(alldata.cancers));
if strcmpi(score_method,'min')
    gind = any(sel_data > min_gene_score, 2);
%     cind = any(sel_data > min_gene_score, 1);  % optional cancer type filtration
elseif strcmpi(score_method,'mean')
    gind = mean(sel_data, 2) > min_gene_score;
elseif strcmpi(score_method,'median')
    gind = median(sel_data, 2) > min_gene_score;
end

% generate heatmap
cmap = flipud(custom_cmap('whitemagma'));
figure
genHeatMap(sel_data(gind,cind), alldata.cancers(cind), alldata.genes(gind), 'both', 'euclid', cmap, [], 'k');

% generate boxplot
[~,sort_gene_ind] = sort(median(sel_data(gind,:),2));
plotdata = sel_data(gind,:)';
plotdata = plotdata(:, sort_gene_ind);
genenames = alldata.genes(gind);
genenames = genenames(sort_gene_ind);
figure
boxplot(plotdata, 'PlotStyle', 'compact', 'Orientation', 'horizontal', 'Colors', 'k');
set(gca,'YTickLabel', genenames, 'YTick', 1:numel(genenames));


%% Compare model ROC AUC scores among different cancer types and mutations

%*********** specify parameters ***********
group_by = 'muts';  % 'cancers' or 'muts'
comb_method = 'mean';  % 'mean' or 'median' (how to combine scores of each model type)
score_types = {'all'}; %{'ExtraTreesClassifier','RandomForestClassifier','XGBClassifier','SVClinear'};
sort_boxes = true;  % 'true' or 'false' (sort boxes in boxplot by median scores)
%******************************************

% set score type if 'all' is selected
if strcmpi(score_types,'all')
    score_inds = true(size(alldata.modelscores.names));
else
    score_inds = ismember(alldata.modelscores.names,score_types);
end

% separate cancer types and mutations in labels
tmp = split(alldata.cancers,'-mut');
cancers = tmp(:,1);
if size(tmp,2) > 1
    muts = tmp(:,2);
else
    muts = [];
end

% determine unique cancers and mutations
uniq_cancers = unique(cancers);
uniq_muts = unique(muts);

% obtain scores for each cancer-mutation combination
if strcmpi(comb_method,'mean')
    mscores = mean(alldata.modelscores.rocauc(score_inds,:),1)';
elseif strcmpi(comb_method,'median')
    mscores = median(alldata.modelscores.rocauc(score_inds,:),1)';
end

% specify grouping parameter
if strcmpi(group_by,'cancers')
    group_data = cancers;
    group_names = uniq_cancers;
elseif strcmpi(group_by,'muts')
    group_data = muts;
    group_names = uniq_muts;
end

% sort groups by median score, if specified
if ( sort_boxes )
    [~,sort_ind] = sort(-cellfun(@(x) median(mscores(ismember(group_data,x))), group_names));
else
    sort_ind = 1:numel(group_names);
end

% generate boxplot
boxplot(mscores,group_data,'GroupOrder',group_names(sort_ind));
set(gca,'XTickLabel',group_names(sort_ind),'XTickLabelRotation',90);

% clear unneeded intermediate variables
clearvars -except alldata










%% Perform Gene Set Analysis

%****** specify data to analyze ******
gene_names = alldata.genes;
cancers = alldata.cancers;
gstats = alldata.Average;
%*************************************

sig_GSnames = {};
sig_GSsizes = [];
sig_padj_nondir = [];
allGOsets = {'H';'C2CP_BIO';'C2CP_KEGG';'C2CP_REAC';'C5BP';'C5CC';'C5MF'};
for k = 1:length(allGOsets)

% load and process gene set list
humanGO = load_human_GO_terms('MSigDB',allGOsets{k});
GSlist = humanGO;
% GSlist = ppi_GS_names(ppi_GS_scores >= 900,:);
% keep_ind = (~ismember(ppi_GS(:,1),psn_genes)) & ismember(ppi_GS(:,2),psn_genes) & (ppi_GS_scores >= 900);
% GSlist = ppi_GS_names(keep_ind,:);

% remove rows without specified genes
GSlist(~ismember(GSlist(:,2),gene_names),:) = [];

% remove repeated entries in GSlist
[~,gsnums] = ismember(GSlist,GSlist);
[~,uniq_ind] = unique(gsnums,'rows');
GSlist = GSlist(uniq_ind,:);

% remove gene sets with only one entry
[gsnames,gsfreq] = cellfreq(GSlist(:,1));
gscut = gsnames(gsfreq == 1);
GSlist(ismember(GSlist(:,1),gscut),:) = [];

% merge similar GSs
[GSlist,mergedGSinfo] = consolidate_gene_sets(GSlist,'union',0);


% run gene set analysis for each cancer type
padj_nondir = []; p_nondir = [];
GSR = {};
for i = 1:length(cancers)
    
    fprintf('\n\nProcessing cancer type %u of %u.\n\n',i,length(cancers));
    
    [GSAres,gene_sets] = GeneSetAnalysis(gene_names,gstats(:,i),[],GSlist,'wilcoxon',10000,[5 Inf],'other');
    
    GSR{i,1} = GSAres;
    
    [~,ind] = ismember({'padj_nondir';'p_nondir'},GSAres(1,:));
    padj_nondir(:,i) = cell2mat(GSAres(2:end,ind(1)));
    p_nondir(:,i) = cell2mat(GSAres(2:end,ind(2)));
    
end
GSnames = GSR{1}(2:end,1);
GSsizes = cell2mat(GSR{1}(2:end,2));
GSRcols = GSR{1}(1,:)';

% filter out nonsignificant gene sets
ind = any(padj_nondir < 0.05,2);
% ind = true(length(GSnames),1);


sig_GSnames = [sig_GSnames;GSnames(ind)];
sig_GSsizes = [sig_GSsizes;GSsizes(ind)];
sig_padj_nondir = [sig_padj_nondir;padj_nondir(ind,:)];
end

% plot results
cg = RNAseq_heatmap(-log10(padj_nondir(ind,:)),cancers,GSnames(ind),'none','integer','hot');


% % identify repeated gene sets, and keep only the last appearance of each
% cut_ind = [];
% for i = 1:length(sig_GSnames)
%     if sum(ismember(sig_GSnames(i:end),sig_GSnames(i))) > 1
%         cut_ind = [cut_ind;i];
%     end
% end
% sig_GSnames(cut_ind) = [];
% sig_GSsizes(cut_ind) = [];
% sig_padj_nondir(cut_ind,:) = [];


%% Group by cancer type or mutation

% option to merge COAD and READ as single cancer type
merge_COAD_READ = false;

% extract list of cancer types and mutations
cancer_muts = cellsplit(alldata.cancers,'mut');
cancers = cancer_muts(:,1);
muts = cancer_muts(:,2);

% merge COAD and READ if requested
if ( merge_COAD_READ )
    cancers(ismember(cancers,{'COAD','READ'})) = {'COADREAD'};
end

% determine frequency of each mutation and cancer
[c_uniq,c_freq] = cellfreq(cancers);
[m_uniq,m_freq] = cellfreq(muts);

% determine mean score and rank for each mutation and cancer
clear merged_cancers merged_muts
merged_cancers.cancers = c_uniq;
merged_muts.muts = m_uniq;
score_fields = fields(alldata);
score_fields(ismember(score_fields,{'genes','cancers','modelscores'})) = [];
for i = 1:length(c_uniq)
    ind = ismember(cancers,c_uniq(i));
    for j = 1:length(score_fields)
        merged_cancers.(score_fields{j})(:,i) = median(alldata.(score_fields{j})(:,ind),2);
    end
end
for i = 1:length(m_uniq)
    ind = ismember(muts,m_uniq(i));
    for j = 1:length(score_fields)
        merged_muts.(score_fields{j})(:,i) = median(alldata.(score_fields{j})(:,ind),2);
    end
end


% assemble cancer-mutation maps for each cancer type
if ~exist('clin_data','var')
    if ~exist('mut_data','var')
        load_mutation_data  % load mutation data if not already loaded
    end
    load_clinical_data  % load clinical data if not already loaded
end

cancer_mut_map = {};  % initialize data structure
for i = 1:length(c_uniq)  % iterate through cancer types
    
    % get list of mutations associated with current cancer type
    m = cancer_muts(ismember(cancer_muts(:,1),c_uniq(i)),2);
    cancer_mut_map.(c_uniq{i}).genes = m;
    
    % determine clinical data indices matching cancer type
    cind = cellfind(c_uniq(i),clin_data.project_id);
    
    % iterate through mutated genes for current cancer type
    add_mut_data = zeros(sum(cind),length(m));
    for j = 1:length(m)
        gind = ismember(clin_data.mut_genes,m(j));  % determine clinical data gene index matching mutated gene
        add_mut_data(:,j) = clin_data.mut_data(cind,gind);  % extract mutation data for current gene and cancer type
    end
    cancer_mut_map.(c_uniq{i}).data = add_mut_data;

end

% visualize the similarity among mutated/non-mutated groups within a cancer type
cancer_type = 'STAD';
mutdist = class_distance(cancer_mut_map.(cancer_type).data');
cg = RNAseq_heatmap(mutdist,cancer_mut_map.(cancer_type).genes,cancer_mut_map.(cancer_type).genes,'none','corr','hot');
pcolor_from_clustergram(cg,mutdist,cancer_mut_map.(cancer_type).genes,cancer_mut_map.(cancer_type).genes,[],'none','jet',[0 1]);

% compare similarity of mutation profile to correlation among gene-level scores
cind = cellfind(strcat(cancer_type,'mut'),alldata.cancers);  % collect scores for only selected cancer type
gind = any(alldata.goodavg(:,cind) > 0.8,2);  % only include genes scoring above some threshold score for the selected cancer type
score_corr = corr(alldata.goodavg(gind,cind));  % calculate the correlation among gene-level scores for all pairs of mutations
if ~all(strcmp(muts(cind),cancer_mut_map.(cancer_type).genes))
    % verify that the mutation names and their ordering match between the score correlation matrix and the mutation distance matrix.
    error('mutation map and gene-level score matrices are not properly aligned.');
end
x = 1 - mutdist(triu(~eye(length(mutdist))));
y = score_corr(triu(~eye(length(score_corr))));
plot(x,y,'o'); xlabel('mutation similarity'); ylabel('correlation of gene-level scores');


%% Extract relevant portion of PPI network

% set parameters
genes = alldata.genes;
map_data = [alldata.Average,mean(alldata.Average,2)];
map_data_labels = [alldata.cancers;{'mean_all_cancers'}];

% load 'InBio' PPI network data
ppi_data_source = 'inbio';
load_ppi_network

% identify subset of PPI network involving map_genes
keep_prots = uniprot2name(ismember(uniprot2name(:,2),genes),1);
ind = any(ismember(ppi,keep_prots),2);  % keep PPIs involving at least one gene in list
% ind = all(ismember(ppi,keep_prots),2);  % keep PPIs involving only genes in list

% extract PPI subset
ppi = ppi(ind,:);
ppi_evidence = ppi_evidence(ind);
ppi_scores = ppi_scores(ind);
ppi_init_scores = ppi_init_scores(ind);

% map data to ppi network nodes (proteins)
ppi_nodes = unique(ppi(:));
ppi_nodes_name = ppi_nodes;
ppi_nodes_data = NaN(length(ppi_nodes),length(map_data_labels));
ppi_nodes_type = repmat({'other'},size(ppi_nodes));
h = waitbar(0,'Processing...');
for i = 1:length(ppi_nodes)
    % map gene names and data to protein nodes
    % NOTE: If multiple genes correspond to a protein (UniProtID), combine
    %       those gene names and their associated data (mean) into a single
    %       node label and score respectively.
    waitbar(i/length(ppi_nodes),h);
    ind = ismember(uniprot2name(:,1),ppi_nodes(i));
    if sum(ind) > 0
        ppi_nodes_name(i,1) = cellcat(uniprot2name(ind,2)','/');
        if any(ismember(uniprot2name(ind,2),alldata.genes))
            ppi_nodes_type(i) = {'PSN_protein'};
        end
        ind = ismember(genes,uniprot2name(ismember(uniprot2name(:,1),ppi_nodes(i)),2));
        ppi_nodes_data(i,:) = mean(map_data(ind,:),1);
    end
end
close(h);
[~,ind] = ismember(ppi,ppi_nodes);
ppi_names = ppi_nodes_name(ind);


% write txt files for cytoscape
network_info = [{'prot1','prot2','ppi_score'};[ppi_names,num2cell(ppi_scores)]];
writecell(network_info,'ppi_network_data.txt',true,'\t','%s\t%s\t%1.3f\n');

node_info = [[{'node_name','node_type'},map_data_labels'];[ppi_nodes_name,ppi_nodes_type,num2cell(ppi_nodes_data)]];
nodeInfoFormatSpec = cellcat([{'%s','%s'},repmat({'%1.8f'},1,size(ppi_nodes_data,2))],'\t');
nodeInfoFormatSpec = [nodeInfoFormatSpec{1},'\n'];
writecell(node_info,'ppi_node_data.txt',true,'\t',nodeInfoFormatSpec);


%% Compare results to text-mined gene list


% select term(s) of interest from below:
%   > bladder_carcinoma
%   > breast_cancer
%   > breast_carcinoma
%   > uterine_cancer
%   > uterine_carcinoma
%   > uterine_endometrial_carcinoma
searchterm = {'breast_carcinoma'};

% import and process ranked gene list for term of interest
if iscell(searchterm)
    g_mined = {};
    for i = 1:length(searchterm)
        tmp = importdata(strcat('/Users/jonrob/Documents/PostDoc/BigDataProject/results/text_mining/beegle_search_',searchterm{i},'.tsv'));
        tmp = cellsplit(tmp,{'"\t"'},true);
        g_mined = union(g_mined,tmp(2:end,2));
    end
else 
    tmp = importdata(strcat('/Users/jonrob/Documents/PostDoc/BigDataProject/results/text_mining/beegle_search_',searchterm,'.tsv'));
    tmp = cellsplit(tmp,{'"\t"'},true);
    g_mined = unique(tmp(2:end,2));
end
clear tmp



% test enrichment of high-scoring genes in text-mined gene list

%****** specify data to analyze ******
% gscores = alldata.goodavg(:,1);
cancer_type = 'BRCA';
%*************************************

GSlist = [repmat({'.'},length(g_mined),1),g_mined];


padj_nondir = []; p_nondir = [];
for i = 1:length(alldata.modeltypes)
    
    fprintf('\n\nProcessing model type %u of %u.\n\n',i,length(alldata.modeltypes));
    
    [GSAres,gene_sets] = GeneSetAnalysis(alldata.genes,alldata.(cancer_type)(:,i),[],GSlist,'wilcoxon',10000,[1 Inf],'other');
    
    GSR{i,1} = GSAres;
    
    [~,ind] = ismember({'padj_nondir';'p_nondir'},GSAres(1,:));
    padj_nondir(:,i) = cell2mat(GSAres(2:end,ind(1)));
    p_nondir(:,i) = cell2mat(GSAres(2:end,ind(2)));
    
end




