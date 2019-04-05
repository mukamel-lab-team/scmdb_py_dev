var trace_3d, layout_3d;
var updates_3d = [];
var groups_3d = [];
var num_colors = 0;

// HTML5 local storage with expiration
// https://gist.github.com/anhang/1096149
var storage = {
	save : function(key, jsonData, expirationMin){
		if (typeof (Storage) === "undefined"){return false;}
        let expirationMS = expirationMin * 60 * 1000;
		let record = {value: JSON.stringify(jsonData), timestamp: new Date().getTime() + expirationMS}
		localStorage.setItem(key, JSON.stringify(record));
		return jsonData;
	},
	load : function(key){
		if (typeof (Storage) === "undefined"){return false;}
		let record = JSON.parse(localStorage.getItem(key));
		if (!record){return false;}
		return (new Date().getTime() < record.timestamp && JSON.parse(record.value));
	}
}

function delayLoad(f) {
    setTimeout(f, 50);
}

function save3DData(trace, layout){
    trace_3d = trace;
    layout_3d = layout;
}

function storeUpdate(update, group, empty=false) {
    if (empty === false){
        updates_3d.push(update);
        groups_3d.push(group);
    }
    else {
        updates_3d = [];
        groups_3d = [];
    }
}

function display3DPlotToggle() {
    if ($('#toggle-3d').prop('checked')){
        $('#loading-3d-plot').html("Loading..");
        Plotly.newPlot("plot-3d-cluster", Object.values(trace_3d), layout_3d);
        if($('#methylation-tsne-grouping option:selected').val() === 'biosample'){
            for(i = 0; i < groups_3d.length; i++){
                Plotly.restyle("plot-3d-cluster", updates_3d[i], groups_3d[i]);
            }
        }
        $('#loading-3d-plot').html("");
        $('#plot-2d-cluster').hide();
        $('#plot-3d-cluster-div').show();
    }
    else {
        Plotly.purge("plot-3d-cluster");
        $('#plot-2d-cluster').show();
        $('#plot-3d-cluster-div').hide();
    }
}

function getMax(arr, prop) {
    let max = 0;
    for (let key in arr) {
        if(parseInt(arr[key][prop]) > max)
            max = arr[key][prop];
    }
    return max;
}

function generateBrowserURL(gene) {
    let base = 'http://brainome.ucsd.edu/annoj_private/CEMBA/cemba.php?ens='+ensemble.replace(/Ens/g, "");

    const chrom = gene.chr.replace(/^\D+/g, "");

    if (gene.strand === '+') {
        var position = gene.start;
    } else {
        var position = gene.end;
    }
    return base+'&assembly='+chrom+'&position='+position;
}

function initGeneNameSearch() {
    geneNameSelector = $('#geneName').select2({
        placeholder: 'Search..',
        allowClear: true,
        ajax: {
            url: './gene/names',
            dataType: 'json',
            delay: 500,
            data: function(params) {
                return {
                    q: params.term
                };
            },
            processResults: function(data) {
                geneSearchCache = data;
                return {
                    results: $.map(data, function(gene) {
                        return {
                            text: gene.gene_name,
                            id: gene.gene_id
                        }
                    })
                }
            },
            cache: true
        },
        minimumInputLength: 1
    });

    //Initialise selector
    let defaultGene = storage.load('lastViewedGenes');
    if (!defaultGene || defaultGene.length === 0) {
        //no entry or browser does not support localStorage, set default to GAD2
        defaultGene = [{gene_name: 'GAD2', gene_id: 'ENSMUSG00000026787'}, {gene_name: 'GAD2', gene_id: 'ENSG00000136750'}];
    }
    console.log('Default gene: '+defaultGene[0])

    if(defaultGene !== []){
        let numGenes = defaultGene.length;
        genes_query = "";
        for (i = 0; i < numGenes; i++) {
            genes_query += defaultGene[i].gene_id + '+';
        }
        genes_query = genes_query.slice(0, -1);
        $.ajax({
            url: './gene/id?q='+genes_query,
            dataType: 'json',
            async: false,
            success: function(data) {
                if (data.length !== 0) {
                    $.each(data, function(i, gene) {
                        let option = new Option(gene.gene_name, gene.gene_id, true, true);
                        geneNameSelector.append(option);
                    });
                    $('#epiBrowserLink').attr('href', generateBrowserURL(data[0]));  // Show the first gene in the list
                }
            }
        });
        updateGeneElements();
    }

    //Hacky way of making sure select2 does not automatically sort tags
    geneNameSelector.on('select2:select', function(evt) {
        var element = evt.params.data.element;
        var $element = $(element);

        $element.detach();
        $(this).append($element);
        $(this).trigger("change");
    });
}

function initGeneModules() {
     geneModuleSelector = $('#geneModulesSelect').select2({
        placeholder: 'Select..',
        allowClear: true,
        minimumResultsForSearch: Infinity
    });

    $.getJSON({
        url: './gene/modules',
        success: function(data){
            data.forEach(function(gene) {
                let option = new Option(gene.module, gene.module, false, false);
                geneModuleSelector.append(option);
            });
        }
    });
}

function updateSearchWithModules(module) {
	$.getJSON({
		url: './gene/modules?q='+module.id,
		success: function (data) {
			data.forEach(function(gene) {
				let option = new Option(gene.gene_name, gene.gene_id, true, true);
				geneNameSelector.append(option);
			});
		}
	});
}


// Drop down options for tSNE (methylation) plot //
function populateMethylationTSNEDropdowns(data) {
    window.global_all_methylation_tsne_settings = data['all_tsne_settings'];
    window.global_all_methylation_clustering_settings = data['all_clustering_settings'];
    window.global_all_methylation_clustering_settings2 = data['all_clustering_settings2'];

    $.each(data["methylation_metadata_fields"], function(key, val) {
        $("#methylation-tsne-grouping").append(
            $('<option></option>').val(val).text(val)
        );
    });
    $.each(data["snATAC_metadata_fields"], function(key, val) {
        $("#snATAC-tsne-grouping").append(
            $('<option></option>').val(val).text(val)
        );
    });
    $.each(data["RNA_metadata_fields"], function(key, val) {
        $("#RNA-tsne-grouping").append(
            $('<option></option>').val(val).text(val)
        );
    });
    
    $.each(data["tsne_dimensions"], function(key, val) {
        $("#methylation-tsne-dimensions").append(
            $('<option></option>').val(val).text(val)
        );
    });
    dynamicMethylationTSNEOptions_dimensions();

    // $.each(data['tsne_methylation'], function(key, val) {
    //     $(".methylation-tsne-methylation").append(
    //         $('<option></option>').val(val).text(val)
    //     );
    // });
    // $.each(data["tsne_perplexity"], function(key, val) {
    //     $(".methylation-tsne-perplexity").append(
    //         $('<option></option>').val(val).text(val)
    //     );
    // });

    $.each(data["clustering_algorithms"], function(key, val) {
        $(".methylation-clustering-algorithms").append(
            $('<option></option>').val(val).text(val)
        );
    });
    dynamicMethylationClusteringOptions_algorithm();
    /*
    $.each(data["clustering_methylation"], function(key, val) {
        $(".methylation-clustering-methylation").append(
            $('<option></option>').val(val).text(val)
        );
    });
    $.each(data["clustering_npc"], function(key, val) {
        $(".methylation-clustering-npc").append(
            $('<option></option>').val(val).text(val)
        );
    });
    $.each(data["clustering_k"], function(key, val) {
        $(".methylation-clustering-k").append(
            $('<option></option>').val(val).text(val)
        );
    });
        */
}

// Options for tSNE (snATAC) plot //
function populatesnATACTSNEDropdowns() {
    $.ajax({
        url: '/snATAC_tsne_options/'+ensemble,
        dataType: 'json',
        async: false,
        success: function(data) {
            window.global_all_snATAC_tsne_settings = data['all_tsne_settings'];
            window.global_all_snATAC_clustering_settings = data['all_clustering_settings'];
            $.each(data["tsne_dimensions"], function(key, val) {
                $(".snATAC-tsne-dimensions").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["tsne_perplexity"], function(key, val) {
                $(".snATAC-tsne-perplexity").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["clustering_algorithms"], function(key, val) {
                $(".snATAC-clustering-algorithms").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["clustering_npc"], function(key, val) {
                $(".snATAC-clustering-npc").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["clustering_k"], function(key, val) {
                $(".snATAC-clustering-k").append(
                    $('<option></option>').val(val).text(val)
                );
            });
        }
    });
}

// Options for tSNE (RNA) plot //
function populateRNATSNEDropdowns() {
    $.ajax({
        url: '/RNA_tsne_options/'+ensemble,
        dataType: 'json',
        async: false,
        success: function(data) {
            window.global_all_RNA_tsne_settings = data['all_tsne_settings'];
            window.global_all_RNA_clustering_settings = data['all_clustering_settings'];
            $.each(data["tsne_dimensions"], function(key, val) {
                $(".RNA-tsne-dimensions").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["tsne_perplexity"], function(key, val) {
                $(".RNA-tsne-perplexity").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["clustering_algorithms"], function(key, val) {
                $(".RNA-clustering-algorithms").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["clustering_npc"], function(key, val) {
                $(".RNA-clustering-npc").append(
                    $('<option></option>').val(val).text(val)
                );
            });
            $.each(data["clustering_k"], function(key, val) {
                $(".RNA-clustering-k").append(
                    $('<option></option>').val(val).text(val)
                );
            });
        }
    });
}

function dynamicMethylationTSNEOptions_dimensions(matching_tsne_options = []) {
    try {  // In case the ensemble does not contain methylation cells. TODO: Fix this so that it uses the tSNE options for ATAC
        if (matching_tsne_options.length === 0) {
            const regex = new RegExp("_ndim" + $(".methylation-tsne-dimensions").val() + "_\\w+");
            for (let i = 0; i < global_all_methylation_tsne_settings.length; i++) {
                if (global_all_methylation_tsne_settings[i].match(regex) !== null) {
                    matching_tsne_options.push(global_all_methylation_tsne_settings[i]);
                }
            }
        }

        let methylation_set = new Set();
        for (let i = 0; i < matching_tsne_options.length; i++) {
            let methylation_type = matching_tsne_options[i].split('_')[0];
            methylation_set.add(methylation_type);
        }

        let methylation_list = [...methylation_set];

        $(".methylation-tsne-methylation").empty();
        $(".methylation-tsne-perplexity").empty();
        for (let i = 0; i < methylation_list.length; i++) {
            if (methylation_list[i] === "mCHmCG") {
                $(".methylation-tsne-methylation").append(
                    $("<option selected></option>").val(methylation_list[i]).text(methylation_list[i])
                );
            }
            else {
                $(".methylation-tsne-methylation").append(
                    $("<option></option>").val(methylation_list[i]).text(methylation_list[i])
                );
            }
        }
        dynamicMethylationTSNEOptions_methylation(matching_tsne_options);
    } catch(error) {}
}

// Updating tsne methylation Levels
function dynamicMethylationTSNEOptions_methylation(matching_tsne_options = []) {
    try { 
        if (matching_tsne_options.length === 0) {
            const regex = new RegExp('^'+$(".methylation-tsne-methylation").val()+"_ndim"+$(".methylation-tsne-dimensions").val()+ "_\\w+");
    		// console.log(regex);
            for (let i = 0; i < global_all_methylation_tsne_settings.length; i++) {
                if (global_all_methylation_tsne_settings[i].match(regex) !== null) {
                    matching_tsne_options.push(global_all_methylation_tsne_settings[i]);
                }
            }
        }

        let perplexity_set = new Set();
        for (let i = 0; i < matching_tsne_options.length; i++) {
            let dimensions = matching_tsne_options[i].split('_')[2];
            perplexity_set.add(dimensions.replace("perp",""));
        }
        let perplexity_list = [...perplexity_set];

        $(".methylation-tsne-perplexity").empty();
        for (let i = 0; i < perplexity_list.length; i++) {
            $(".methylation-tsne-perplexity").append(
                $("<option></option>").val(perplexity_list[i]).text(perplexity_list[i])
            );
        }
    } catch(error) {}
}

function dynamicMethylationClusteringOptions_algorithm() {
    try {
        var matching_clustering_options = [];
        const regex = new RegExp('_' + $(".methylation-clustering-algorithms").val() + "_\\w+");
    
        for (let i = 0; i < global_all_methylation_clustering_settings.length; i++) {
            if (global_all_methylation_clustering_settings[i].match(regex) !== null) {
                matching_clustering_options.push(global_all_methylation_clustering_settings[i]);
            }
        }
    
        let methylation_set = new Set();
        for (let i = 0; i < matching_clustering_options.length; i++) {
            let methylation_type = matching_clustering_options[i].split('_')[0];
            methylation_set.add(methylation_type);
        }
        let methylation_list = [...methylation_set];
    
        $(".methylation-clustering-methylation").empty();
        for (let i = 0; i < methylation_list.length; i++) {
            if (methylation_list[i] === 'mCHmCG') {
                $(".methylation-clustering-methylation").append(
                    $("<option selected></option>").val(methylation_list[i]).text(methylation_list[i])
                );
            }
            else {
                $(".methylation-clustering-methylation").append(
                    $("<option></option>").val(methylation_list[i]).text(methylation_list[i])
                );
            }
        }
    
        dynamicMethylationClusteringOptions_methylation();
    } catch(error) {}
}

// No longer necessary
function dynamicMethylationClusteringOptions_methylation(matching_clustering_options = []) {
    try {
        if (matching_clustering_options.length === 0) {
            var matching_clustering_options = [];
            const regex = new RegExp('^'+$(".methylation-clustering-methylation").val()+'_'+$(".methylation-clustering-algorithms").val() + "_\\w+");
    
            for (let i = 0; i < global_all_methylation_clustering_settings.length; i++) {
                if (global_all_methylation_clustering_settings[i].match(regex) !== null) {
                    matching_clustering_options.push(global_all_methylation_clustering_settings[i]);
                }
            }
        }
    
        let npc_set = new Set();
        for (let i = 0; i < matching_clustering_options.length; i++) {
            let npc = matching_clustering_options[i].split('_')[2].replace("npc","");
            npc_set.add(npc);
        }
        let npc_list = [...npc_set];
    
        $(".methylation-clustering-npc").empty();
        for (let i = 0; i < npc_list.length; i++) {
            $(".methylation-clustering-npc").append(
                $("<option></option>").val(npc_list[i]).text(npc_list[i])
            );
        }
    
        dynamicMethylationClusteringOptions_npc(matching_clustering_options);
    } catch(error) {}
}


function dynamicMethylationClusteringOptions_npc(matching_clustering_options = []) {
    if (matching_clustering_options.length === 0) {
        var matching_clustering_options = [];
        const regex = new RegExp('^'+$(".methylation-clustering-methylation").val()+'_'+$(".methylation-clustering-algorithms").val()+'_npc50' + "_\\w+");

        for (let i = 0; i < global_all_methylation_clustering_settings.length; i++) {
            if (global_all_methylation_clustering_settings[i].match(regex) !== null) {
                matching_clustering_options.push(global_all_methylation_clustering_settings[i]);
                //matching_clustering_options.push(i);
            }
        }
    }

    k_list = [];
    k_clusters_list = [];
    for (let i = 0; i < matching_clustering_options.length; i++) {
        let k = matching_clustering_options[i].split('_')[3].replace("k","");
        let clusters = global_all_methylation_clustering_settings2[matching_clustering_options[i]];
        k_list.push(k);
        k_clusters_list.push(k+' ('+clusters+' clusters)');
    }

    $(".methylation-clustering-k").empty();
    for (let i = 0; i < k_list.length; i++) {
        if (k_list[i] === "20") {
            $(".methylation-clustering-k").append(
                $("<option selected></option>").val(k_list[i]).text(k_clusters_list[i])
            );
        }
        else {
            $(".methylation-clustering-k").append(
                $("<option></option>").val(k_list[i]).text(k_clusters_list[i])
            );
        }
    }
}

function dynamicsnATACTSNEOptions_dimensions(matching_tsne_options = []) {
    if (matching_tsne_options.length === 0) {
        var matching_tsne_options = [];
        const regex = new RegExp("ATAC_ndim" + $(".snATAC-tsne-dimensions").val() + "_\\w+");
        for (let i = 0; i < global_all_snATAC_tsne_settings.length; i++) {
            if (global_all_snATAC_tsne_settings[i].match(regex) !== null) {
                matching_tsne_options.push(global_all_snATAC_tsne_settings[i]);
            }
        }
    }

    let perplexity_set = new Set();
    for (let i = 0; i < matching_tsne_options.length; i++) {
        let perplexity = matching_tsne_options[i].split('_')[2];
        perplexity_set.add(perplexity.replace("perp",""));
    }
    let perplexity_list = [...perplexity_set];

    $(".snATAC-tsne-perplexity").empty();
    for (let i = 0; i < perplexity_list.length; i++) {
        $(".snATAC-tsne-perplexity").append(
            $("<option></option>").val(perplexity_list[i]).text(perplexity_list[i])
        );
    }
}


function dynamicsnATACClusteringOptions_algorithm(matching_clustering_options = []) {
    if (matching_clustering_options.length === 0) {
        var matching_clustering_options = [];
        const regex = new RegExp('ATAC_'+$(".methylation-clustering-algorithms").val() + "_\\w+");

        for (let i = 0; i < global_all_snATAC_clustering_settings.length; i++) {
            if (global_all_snATAC_clustering_settings[i].match(regex) !== null) {
                matching_clustering_options.push(global_all_snATAC_clustering_settings[i]);
            }
        }
    }

    let npc_set = new Set();
    for (let i = 0; i < matching_clustering_options.length; i++) {
        let npc = matching_clustering_options[i].split('_')[2].replace("npc","");
        npc_set.add(npc);
    }
    let npc_list = [...npc_set];

    $(".snATAC-clustering-npc").empty();
    for (let i = 0; i < npc_list.length; i++) {
        $(".snATAC-clustering-npc").append(
            $("<option></option>").val(npc_list[i]).text(npc_list[i])
        );
    }

    dynamicsnATACClusteringOptions_npc(matching_clustering_options);
}

function dynamicsnATACClusteringOptions_npc(matching_clustering_options = []) {
    if (matching_clustering_options.length === 0) {
        var matching_clustering_options = [];
        const regex = new RegExp('ATAC_'+$(".snATAC-clustering-algorithms").val()+'_npc'+$(".snATAC-clustering-npc").val() + "_\\w+");

        for (let i = 0; i < global_all_snATAC_clustering_settings.length; i++) {
            if (global_all_snATAC_clustering_settings[i].match(regex) !== null) {
                matching_clustering_options.push(global_all_snATAC_clustering_settings[i]);
            }
        }
    }

    let k_set = new Set();
    for (let i = 0; i < matching_clustering_options.length; i++) {
        let k = matching_clustering_options[i].split('_')[3].replace("k","");
        k_set.add(k);
    }
    let k_list = [...k_set];

    $(".snATAC-clustering-k").empty();
    for (let i = 0; i < k_list.length; i++) {
        $(".snATAC-clustering-k").append(
            $("<option></option>").val(k_list[i]).text(k_list[i])
        );
    }
}

function updateGeneElements(updateMCHScatter=true) {
    buttons = document.getElementsByClassName('modebar-btn');
    let geneSelected = $('#geneName option:selected').val();
    if (geneSelected != 'Select..' && $("#geneName").select2('data').length > 0) {
        $('#orthologsToggle').bootstrapToggle('off');
        $('#orthologsToggle').bootstrapToggle('enable');

        let lastViewedGenes = [];
        for(i=0; i<$('#geneName').select2('data').length; i++){
            lastViewedGenes.push({gene_name: $('#geneName option:selected')[i].text, gene_id: $('#geneName option:selected')[i].value});
        }
        if (typeof(Storage) !== 'undefined') {
            storage.save('lastViewedGenes', lastViewedGenes, 30);  // store last viewed genes for 30 minutes
        }
        if (updateMCHScatter){
            if (methylation_data_available === 1) {
                updateMCHScatterPlot();
            }
            if (snATAC_data_available === 1) {
                updatesnATACScatterPlot();
            }
            if (RNA_data_available === 1) {
                updateRNAScatterPlot();
            }
        }
        updateClustersBarPlot();
        $('#epiBrowserLink').removeClass('disabled');
        if($("#geneName").select2('data').length > 1) {
            $('#normalize-heatmap').show();
            $('#methylation-box-heat-normalize-toggle-div').show();
            $('#methylation-box-heat-normalize-toggle').prop('disabled', false);
            updateMethylationHeatmap();
            updateCorrelatingGeneDataTable("");
            // $('#epiBrowserLink').addClass('disabled');
            $("#methylation-box-and-heat").removeClass('col-md-8').addClass('col-md-12');
            $("#methylation-correlated-genes").hide();
            if (RNA_data_available === 1) {
                updateRNAHeatmap();
                $('#RNA-box-heat-normalize-toggle').prop('disabled', false);
            }
        }
        else{
            $('#normalize-heatmap').hide();
            $('#methylation-box-heat-normalize-toggle-div').hide();
            $('#methylation-box-heat-normalize-toggle').prop('disabled', true);
            //updateOrthologToggle();
            updateMCHBoxPlot();
            updateCorrelatingGeneDataTable($('#geneName option:selected').val());
            $("#methylation-box-and-heat").removeClass('col-md-12').addClass('col-md-8');
            $("#methylation-correlated-genes").show();
            if (RNA_data_available === 1) {
                updateRNABoxPlot();
                $('#RNA-box-heat-normalize-toggle').prop('disabled', true);
            }
        }
        $.ajax({
            url: './gene/id?q='+geneSelected,
            dataType: 'json',
            success: function(data) {
                if (typeof(data.gene_name) !== 'undefined' && typeof(data.gene_id) !== 'undefined') {
                    $('#epiBrowserLink').attr('href', generateBrowserURL(data));
                    $('#epiBrowserLink').removeClass('disabled');
                }
            }
        });
    }
}

function updateMCHScatterPlot(onlyUpdatetSNEandClustering=false) {
	// Methylation type
	let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
    let max_points = $('#max-points').val();
    let methylationType = $("#mType").val();
    let methylation_color_percentile_Values = methylation_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let grouping = $('#methylation-tsne-grouping').val();
    let genes_query = "";

	// Used variables: methylation type (tSNE), tSNE perplexity, methylation-tsne-dimensions
    let tsne_setting = $("#methylation-tsne-methylation").val() + "_ndim" + $("#methylation-tsne-dimensions").val() + "_perp" + $("#methylation-tsne-perplexity").val();

	// Used variables: methylationType (clustering), methylation-algorithms(???), k-value
	let clustering = $("#methylation-clustering-methylation").val()+"_"+$("#methylation-clustering-algorithms").val()+"_npc50"+"_k"+$("#methylation-clustering-k").val();
	console.log(clustering);
	// Check the outlier option
    if ($('#methylation_tsneOutlierToggle').prop('checked')) {
        var tsneOutlierOption = 'false';
    } else {
        var tsneOutlierOption = 'true';
    }

	// Update genes if needed
    if (!onlyUpdatetSNEandClustering) {
        for (i = 0; i < genes.length; i++) {
            genes_query += (genes[i].id + "+");
        }
    }
    else {
        let lastViewedGenes = storage.load('lastViewedGenes');
        for (i = 0; i < lastViewedGenes.length; i++) {
            genes_query += (lastViewedGenes[i].gene_id + "+");
        }
    }
    genes_query = genes_query.slice(0,-1);

	// Load the data
    if ($('#geneName option:selected').val() != 'Select..') {
        $.ajax({
        //$.getJSON({
            type: "GET",
            url: './plot/methylation/scatter/'+ensemble+'/'+tsne_setting+'/' +methylationType+ '/'+levelType+'/'+grouping+'/'+clustering+'/'+methylation_color_percentile_Values[0]+'/'+methylation_color_percentile_Values[1]+'/'+tsneOutlierOption+'/'+max_points+'?q='+genes_query,
            beforeSend: function() {
                $("#mch-scatter-loader").show();
                $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr('disabled', true);
                //$("#plot-mch-scatter").html("");
            },
            complete: function() {
                $("#mch-scatter-loader").hide();
            },
            success: function(data) {
                //Plotly.newPlot('plot-mch-scatter', data);
                $('#plot-mch-scatter').html(data);
                $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr('disabled', false);
            }
        });
    }

    $("#methylation-tsne-heading-num-dimensions").text($("#methylation-tsne-dimensions").val() + "D ");
    $("#methylation-tsne-options-heading").text("Methylation: " + $("#methylation-tsne-methylation").val() + ", Perplexity: " + $("#methylation-tsne-perplexity").val());

    $("#methylation-clustering-options-heading").text("Algorithm: " + $("#methylation-clustering-algorithms").val() + ", Methylation: " + $("#methylation-clustering-methylation").val() + ", K-value: " + $("#methylation-clustering-k").val());
}

function updatesnATACScatterPlot(onlyUpdatetSNEandClustering=false) {
    let tsne_settings = "ATAC_ndim"+$("#snATAC-tsne-dimensions").val()+"_perp"+$("#snATAC-tsne-perplexity").val();
    let max_points = $('#max-points').val();
    let grouping = $("#methylation-tsne-grouping").val();
    let clustering_settings = "ATAC_"+$("#snATAC-clustering-algorithms").val()+"_npc"+$("#snATAC-clustering-npc").val()+"_k"+$("#snATAC-clustering-k").val();
    let snATAC_color_percentile_Values = snATAC_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let genes_query = "";
    let smoothing = $('#snATAC-smoothing-toggle').prop('checked');

    if ($('#snATAC_tsneOutlierToggle').prop('checked')) {
        var tsneOutlierOption = 'false';
    } else {
        var tsneOutlierOption = 'true';
    }

    if (!onlyUpdatetSNEandClustering) {
        for (i = 0; i < genes.length; i++) {
            genes_query += (genes[i].id + "+");
        }
    }
    else {
        let lastViewedGenes = storage.load("lastViewedGenes");
        for (i = 0; i < lastViewedGenes.length; i++) {
            genes_query += (lastViewedGenes[i].gene_id + "+");
        }
    }
    genes_query = genes_query.slice(0,-1);

    if ($('#geneName option:selected').val() != 'Select..') {
        $.ajax({
        //$.getJSON({
            type: "GET",
            url: './plot/snATAC/scatter/'+ensemble+'/'+grouping+'/'+snATAC_color_percentile_Values[0]+'/'+snATAC_color_percentile_Values[1]+'/'+tsneOutlierOption+'/'+smoothing+'/'+max_points+'?q='+genes_query,
            beforeSend: function() {
                $("#snATAC-scatter-loader").show();
                $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", true);
            },
            complete: function() {
                $("#snATAC-scatter-loader").hide();
            },
            success: function(data) {
                //Plotly.newPlot('plot-mch-scatter', data);
                $('#plot-snATAC-scatter').html(data);
                $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", false);
            }
        });
    }

    $("#snATAC-tsne-heading-num-dimensions").text($("#snATAC-tsne-dimensions").val() + "D ");
    $("#snATAC-tsne-options-heading").text("Perplexity: " + $("#snATAC-tsne-perplexity").val());
    $("#snATAC-clustering-options-heading").text("Algorithm: " + $("#snATAC-clustering-algorithms").val() + ", # of PCs: " + $("#snATAC-clustering-npc").val() + ", K-value: " + $("#snATAC-clustering-k").val());
}


function updateRNAScatterPlot(onlyUpdatetSNEandClustering=false) {
    let tsne_settings = "ATAC_ndim"+$("#RNA-tsne-dimensions").val()+"_perp"+$("#RNA-tsne-perplexity").val();
    let max_points = $('#max-points').val();
    let grouping = $("#RNA-tsne-grouping").val();
    let clustering_settings = "ATAC_"+$("#RNA-clustering-algorithms").val()+"_npc"+$("#RNA-clustering-npc").val()+"_k"+$("#RNA-clustering-k").val();
    let RNA_color_percentile_Values = RNA_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let genes_query = "";

    if ($('#RNA_tsneOutlierToggle').prop('checked')) {
        var tsneOutlierOption = 'false';
    } else {
        var tsneOutlierOption = 'true';
    }

    if (!onlyUpdatetSNEandClustering) {
        for (i = 0; i < genes.length; i++) {
            genes_query += (genes[i].id + "+");
        }
    }
    else {
        let lastViewedGenes = storage.load("lastViewedGenes");
        for (i = 0; i < lastViewedGenes.length; i++) {
            genes_query += (lastViewedGenes[i].gene_id + "+");
        }
    }
    genes_query = genes_query.slice(0,-1);

    if ($('#geneName option:selected').val() != 'Select..') {
        $.ajax({
        //$.getJSON({
            type: "GET",
            url: './plot/RNA/scatter/'+ensemble+'/'+grouping+'/'+RNA_color_percentile_Values[0]+'/'+RNA_color_percentile_Values[1]+'/'+tsneOutlierOption+'/'+max_points+'?q='+genes_query,
            beforeSend: function() {
                $("#RNA-scatter-loader").show();
                $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", true);
            },
            complete: function() {
                $("#RNA-scatter-loader").hide();
            },
            success: function(data) {
                //Plotly.newPlot('plot-mch-scatter', data);
                $('#plot-RNA-scatter').html(data);
                $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", false);
            }
        });
    }

    $("#RNA-tsne-heading-num-dimensions").text($("#RNA-tsne-dimensions").val() + "D ");
    $("#RNA-tsne-options-heading").text("Perplexity: " + $("#RNA-tsne-perplexity").val());
    $("#RNA-clustering-options-heading").text("Algorithm: " + $("#RNA-clustering-algorithms").val() + ", # of PCs: " + $("#RNA-clustering-npc").val() + ", K-value: " + $("#RNA-clustering-k").val());
}

/*
function updateOrthologToggle() {
    let geneSelected = $('#geneName option:selected').val();
    $.ajax({
        type: "GET",
        url: './gene/orthologs/'+ensemble+'/'+geneSelected,
        success: function(data) {
            if (data.mmu_gID === "" || data.hsa_gID === "") {
                $('#orthologsToggle').bootstrapToggle('off');
                $('#orthologsToggle').bootstrapToggle('disable');
            } else {
                mmu_gID = data.mmu_gID;
                hsa_gID = data.hsa_gID;
                $('#orthologsToggle').bootstrapToggle('enable');
                if ($('#orthologsToggle').prop('checked')) {
                    return updateMCHCombinedBoxPlot(mmu_gID, hsa_gID);
                }
            }
        }
    });
}
*/

function initDataTableClick() {
    $('#corrGeneTable tbody').on('click', 'tr', function () {
        let id = $(this).attr('id');
        $.getJSON({
            url: './gene/id?q='+id,
            success: function (data) {
                let option = new Option(data[0].gene_name, data[0].gene_id, true, true);
                for(let i=0; i < $("#geneName").select2('data').length; i++) {
                    if($("#geneName").select2('data')[i].id === option.value) {
                        return;
                    }
                }
                geneNameSelector.val(null).trigger("change.select2"); // Clear gene search bar
                geneNameSelector.append(option);
                $('#epiBrowserLink').attr('href', generateBrowserURL(data[0]));
                $('#epiBrowserLink').removeClass('disabled');
                updateGeneElements();
            }
        });
    });
}

function updateCorrelatingGeneDataTable(geneSelected) {
    if (geneSelected !== 'Select..' && geneSelected !== "") {
        table = $('#corrGeneTable').DataTable( {
            "destroy": true,
            "ordering": false,
            "lengthChange": false,
            "dom": "<'col-sm-12'<f>>" +
                    "<<t>>" +
                    "<'col-sm-12'<i>>" +
                    "<'col-sm-12'<p>>",
            "pagingType": "simple",
            "ajax": {
                "url": "./gene/corr/"+ensemble+"/"+geneSelected,
                "dataSrc": "",
            },
            "rowId": 'gene_id',
            "columns": [
                { "data": "rank" },
                { "data": "gene_name" },
                { "data": "correlation" },
            ],
        });
    }
    else {
        table.clear();
    }
}


function initClusterSpecificMarkerGeneTable() {

    let clustering = $("#methylation-clustering-methylation").val()+"_"+$("#methylation-clustering-algorithms").val()+"_npc50"+"_k"+$("#methylation-clustering-k").val();

    let htmlTable = '<table id="clusterMarkerGeneTable" class="display nowrap compact" width="100%"><tbody></tbody></table>';

    // Destroy existing table and create new table.
    if ($('#clusterMarkerGeneTable').length > 0) {
        //$('#clusterMarkerGeneTable').DataTable().destroy(true);
        $('#clusterMarkerGeneTableDiv').html(htmlTable);
    }
    $.ajax({
        type: "GET",
        url: './cluster/marker_genes/'+ensemble+'/'+clustering,
        success: function(data) {
            let columns = [ { mData: "rank", sTitle: "rank" } ];
            $.each(data.columns, function(i, cluster) {
                columns.push({ mData: cluster, sTitle: cluster,
                "defaultContent": "-", });
            });
            clusterMarkerGeneTable = $('#clusterMarkerGeneTable').DataTable( {
                "pageLength": 50,
                "destroy": true,
                "ordering": false,
                "scrollX": "100%",
                "select": true,
                "lengthChange": true,
                "pageLength" : 50,
                "columns": columns,
                "data": data.rows,
            });
            clusterMarkerGeneTable.select.items('column');
            delayLoad(clusterMarkerGeneTable.draw());

            // Fills gene search bar with genes in column when user clicks.
            clusterMarkerGeneTable.on('select', function(e, dt, type, indexes) {
                var data = clusterMarkerGeneTable.columns( indexes ).data();
                genes_query = data[0].join('+');

                $.ajax({
                    type: "GET",
                    url: './gene/names/exact?q='+genes_query,
                    success: function(data) {
                        if (data.length !== 0) {
                            geneNameSelector.val(null).trigger("change.select2"); //Clear gene search bar
                            $.each(data, function(i, gene) {
                                console.log(gene.gene_name);
                                let option = new Option(gene.gene_name, gene.gene_id, true, true);
                                geneNameSelector.append(option);
                            });
                        }
                    }
                });
            });
        }
    });
}

function updateMCHBoxPlot() {
    let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
    let methylationType = $("#mType").val();
    let max_points = $('#max-points').val();
    let geneSelected = $('#geneName option:selected').val();
    let grouping = $('#methylation-tsne-grouping').val();
    let clustering = $("#methylation-clustering-methylation").val()+"_"+$("#methylation-clustering-algorithms").val()+"_npc50"+"_k"+$("#methylation-clustering-k").val();
    if ($('#orthologsToggle').prop('checked')) {
        return updateMCHCombinedBoxPlot(mmu_gID, hsa_gID);
    }
    if ($('#methylation-box-heat-outlierToggle').prop('checked')) {
        var outlierOption = 'outliers';
    } else {
        var outlierOption = 'false';
    }

    $.ajax({
        type: "GET",
        url: './plot/methylation/box/'+ensemble+'/'+methylationType+'/'+geneSelected+'/'+grouping+'/'+clustering+'/'+levelType+'/'+outlierOption+'/'+max_points,
        beforeSend: function() {
            $("#mch-box-loader").show();
            $("#plot-mch-heat").html("");
            $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", true);
        },
        complete: function() {
            $('#mch-box-loader').hide();
        },
        success: function(data) {
            $("#plot-mch-box").html(data);
            $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", false);
        }
    });

}

function updateClustersBarPlot() {
    let grouping = $('#methylation-tsne-grouping').val();
    let normalize = $('#clusters-bar-normalize-toggle').prop('checked');
    let clustering = $("#methylation-clustering-methylation").val()+"_"+$("#methylation-clustering-algorithms").val()+"_npc50"+"_k"+$("#methylation-clustering-k").val();

    $.ajax({
        type: "GET",
        url: './plot/clusters/bar/'+ensemble+'/'+grouping+'/'+clustering+'/'+normalize,
        beforeSend: function() {
            $("#clusters-bar-loader").show();
            $("#plot-clusters-bar").html("");
            $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", true);
        },
        complete: function() {
            $('#clusters-bar-loader').hide();
        },
        success: function(data) {
            $("#plot-clusters-bar").html(data);
            $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", false);
        }
    });
}

function updatesnATACBoxPlot() {
    let geneSelected = $('#geneName option:selected').val();
    let grouping = $('#methylation-tsne-grouping').val();

    if ($('#methylation-box-heat-outlierToggle').prop('checked')) {
        var outlierOption = 'outliers';
    } else {
        var outlierOption = 'false';
    }

    $.ajax({
        type: "GET",
        url: './plot/snATAC/box/'+ensemble+'/'+geneSelected+'/'+grouping+'/'+outlierOption,
        beforeSend: function() {
            // $("#snATAC-box-heat-UpdateBtn").attr("disabled", true);
            $("#snATAC-box-loader").show();
            $("#plot-snATAC-heat").html("");
        },
        complete: function() {
            $("#snATAC-box-loader").hide();
        },
        success: function(data) {
            $('#plot-snATAC-box').html(data);
            // $("#snATAC-box-heat-UpdateBtn").attr("disabled", false);
        }
    });

}

function updateRNABoxPlot() {
    let geneSelected = $('#geneName option:selected').val();
    let grouping = $('#RNA-box-heat-grouping').val();

    if ($('#RNA-box-heat-outlierToggle').prop('checked')) {
        var outlierOption = 'outliers';
    } else {
        var outlierOption = 'false';
    }

    $.ajax({
        type: "GET",
        url: './plot/RNA/box/'+ensemble+'/'+geneSelected+'/'+grouping+'/'+outlierOption,
        beforeSend: function() {
            $("#RNA-box-heat-UpdateBtn").attr("disabled", true);
            $("#RNA-box-loader").show();
            $("#plot-RNA-heat").html("");
        },
        complete: function() {
            $("#RNA-box-loader").hide();
        },
        success: function(data) {
            $('#plot-RNA-box').html(data);
            $("#RNA-box-heat-UpdateBtn").attr("disabled", false);
        }
    });

}

function updateMCHCombinedBoxPlot(mmu_gid, hsa_gid) {
    let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
   let methylationType = $("#mType").val();
    if ($('#methylation-box-heat-outlierToggle').prop('checked')) {
        var outlierOption = 'outliers';
    } else {
        var outlierOption = 'false';
    }

    $.ajax({
        type: "GET",
        url: './plot/box_combined/'+methylationType+'/'+mmu_gid+'/'+hsa_gid+'/'+levelType+'/'+outlierOption,
        success: function(data) {
            $('#plot-mch-heat').html("");
            $('#mch_box_div').addClass("col-md-9");
            $('#gene_table_div').show();
            $('#plot-mch-box').html(data);
        }
    });

}

function updateMethylationHeatmap() {
    let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
    let methylationType = $("#mType").val();
    let methylation_box_color_percentile_Values = methylation_box_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let genes_query = "";
    let grouping = $("#methylation-tsne-grouping").val();
    let clustering = $("#methylation-clustering-methylation").val()+"_"+$("#methylation-clustering-algorithms").val()+"_npc50"+"_k"+$("#methylation-clustering-k").val();

    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    if ($('#methylation-box-heat-normalize-toggle').prop('checked')) {
        var normalize = 'true';
    }
    else {
        var normalize = 'false';
    }
    genes_query = genes_query.slice(0,-1);

    $.ajax({
        type: "GET",
        url: './plot/methylation/heat/'+ensemble+'/'+methylationType+'/'+grouping+'/'+clustering+'/'+levelType+'/'+methylation_box_color_percentile_Values[0]+'/'+methylation_box_color_percentile_Values[1]+'?q='+genes_query+'&normalize='+normalize,
        beforeSend: function() {
            $("#mch-box-loader").show();
            $("#plot-mch-box").html("");
            $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", true);
        },
        complete: function() {
            $("#mch-box-loader").hide();
        },
        success: function(data) {
            $('#gene_table_div').hide();
            $('#mch_box_div').removeClass("col-md-9");
            $('#plot-mch-heat').html(data);
            $("#methylation-tsneUpdateBtn, #methylation-tsneUpdateBtn-top").attr("disabled", false);
            $('#methylation-box-heat-outlierToggle').bootstrapToggle('disable');
        }
    });
}

function updatesnATACHeatmap() {
    let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
    let snATAC_color_percentile_Values = methylation_box_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let genes_query = "";
    let grouping = $("#methylation-tsne-grouping").val();

    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    if ($('#methylation-box-heat-normalize-toggle').prop('checked')) {
        var normalize = 'true';
    }
    else {
        var normalize = 'false';
    }
    genes_query = genes_query.slice(0,-1);

    $.ajax({
        type: "GET",
        url: './plot/snATAC/heat/'+ensemble+'/'+grouping+'/'+snATAC_color_percentile_Values[0]+'/'+snATAC_color_percentile_Values[1]+'?q='+genes_query+'&normalize='+normalize,
        beforeSend: function() {
            $("#snATAC-box-loader").show();
            $("#plot-snATAC-box").html("");
            // $("#snATAC-box-heat-UpdateBtn").attr("disabled", true);
        },
        complete: function() {
            $("#snATAC-box-loader").hide();
        },
        success: function(data) {
            $('#plot-snATAC-heat').html(data);
            $('#methylation-box-heat-outlierToggle').bootstrapToggle('disable');
            // $("#snATAC-box-heat-UpdateBtn").attr("disabled", false);
        }
    });
}

function updateRNAHeatmap() {
    let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
    let RNA_color_percentile_Values = methylation_box_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let genes_query = "";
    let grouping = $("#RNA-box-heat-grouping").val();

    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }

    if ($('#methylation-box-heat-normalize-toggle').prop('checked')) {
        var normalize = 'true';
    }
    else {
        var normalize = 'false';
    }
    genes_query = genes_query.slice(0,-1);

    $.ajax({
        type: "GET",
        url: './plot/RNA/heat/'+ensemble+'/'+grouping+'/'+RNA_color_percentile_Values[0]+'/'+RNA_color_percentile_Values[1]+'?q='+genes_query+'&normalize='+normalize,
        beforeSend: function() {
            $("#RNA-box-loader").show();
            $("#plot-RNA-box").html("");
            $("#RNA-box-heat-UpdateBtn").attr("disabled", true);
        },
        complete: function() {
            $("#RNA-box-loader").hide();
        },
        success: function(data) {
            $('#plot-RNA-heat').html(data);
            $('#RNA-box-heat-outlierToggle').bootstrapToggle('disable');
            $("#RNA-box-heat-UpdateBtn").attr("disabled", false);
        }
    });
}

function updateMethylationHeatmapTwoSpecies() {
    let levelType = $('#methylation-levels').val();
    // let levelType = $('#level').val();
    let methylationType = $("#mType").val();
    let methylation_color_percentile_Values = methylation_color_percentile_Slider.getValue();
    let genes = $("#geneName").select2('data');
    let genes_query = "";

    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    if ($('#methylation-box-heat-normalize-toggle').prop('checked')) {
        var normalize = 'true';
    }
    else {
        var normalize = 'false';
    }
    genes_query = genes_query.slice(0,-1);

    $.ajax({
        type: "GET",
        url: './plot/heat_two_ensemble/'+ensemble+'/'+methylationType+'/'+levelType+'/'+methylation_color_percentile_Values[0]+'/'+methylation_color_percentile_Values[1]+'?q='+genes_query+'&normalize='+normalize,
        success: function(data) {
            $('#plot-mch-box').html("");
            $('#gene_table_div').hide();
            $('#mch_box_div').removeClass("col-md-9");
            $('#plot-mch-heat').html(data);
            $('#methylation-box-heat-outlierToggle').bootstrapToggle('disable');
        }
    });
}
