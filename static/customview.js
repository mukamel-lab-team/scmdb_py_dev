var trace_3d, layout_3d;
var updates_3d = [];
var groups_3d = [];
var num_colors = 0;

// HTML5 local storage with expiration
// https://gist.github.com/anhang/1096149
var storage = {
	save : function(key, jsonData, expirationMin){
		if (typeof (Storage) === "undefined"){return false;}
        var expirationMS = expirationMin * 60 * 1000;
		var record = {value: JSON.stringify(jsonData), timestamp: new Date().getTime() + expirationMS}
		localStorage.setItem(key, JSON.stringify(record));
		return jsonData;
	},
	load : function(key){
		if (typeof (Storage) === "undefined"){return false;}
		var record = JSON.parse(localStorage.getItem(key));
		if (!record){return false;}
		return (new Date().getTime() < record.timestamp && JSON.parse(record.value));
	} }

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
        if($('#tsneGrouping option:selected').val() === 'biosample'){
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
    var max = 0;
    for (var key in arr) {
        if(parseInt(arr[key][prop]) > max)
            max = arr[key][prop];
    }
    return max;
}

function generateBrowserURL(gene) {
    if (ensemble === 'mmu') {
        var base = 'http://brainome.ucsd.edu/annoj/CEMBA/index_mm.html'; // Mouse
    } else {
        var base = 'http://brainome.ucsd.edu/annoj/sc_wgbs/index_hs.html'; // Human
    }

    if (gene.strand === '+') {
        var position = gene.start;
    } else {
        var position = gene.end;
    }
    return base+'?assembly='+gene.chrom+'&position='+position;
}

function initGeneNameSearch() {
    geneNameSelector = $('#geneName').select2({
        placeholder: 'Search..',
        allowClear: true,
        ajax: {
            url: './gene/names/'+ensemble,
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
    var defaultGene = storage.load('lastViewedGenes');
    if (!defaultGene || defaultGene.length === 0) {
        //no entry or browser does not support localStorage, set default to GAD2
        defaultGene = [{gene_name: 'GAD2', gene_id: 'ENSMUSG00000026787'}];
    }
    
    if(defaultGene !== []){
        var numGenes = defaultGene.length;
        for (i = 0; i < numGenes; i++) {
            $.ajax({
                url: './gene/id/'+ensemble+'?q='+defaultGene[i].gene_id,
                dataType: 'json',
                async: false,
                success: function(data) {
                    if (typeof(data.gene_name) !== 'undefined' && typeof(data.gene_id) !== 'undefined') {
                        var option = new Option(data.gene_name, data.gene_id, true, true);
                        geneNameSelector.append(option);
                        if (numGenes === 1) {
                            $('#epiBrowserLink').attr('href', generateBrowserURL(data));
                            $('#epiBrowserLink').removeClass('disabled');
                        }
                    }
                }
            });
        }
        updateGeneElements();
    }
}

function initGeneModules() {
     geneModuleSelector = $('#geneModulesSelect').select2({
        placeholder: 'Select..',
        allowClear: true,
        minimumResultsForSearch: Infinity
    });

    $.getJSON({
        url: './gene/modules/'+ensemble,
        success: function(data){
            data.forEach(function(gene) {
                var option = new Option(gene.module, gene.module, false, false);
                geneModuleSelector.append(option);
            });
        }
    });
}

function updateSearchWithModules(module) {
	$.getJSON({
		url: './gene/modules/'+ensemble+'?q='+module.id,
		success: function (data) {
			data.forEach(function(gene) {
				var option = new Option(gene.gene_name, gene.gene_id, true, true);
				geneNameSelector.append(option);
			});
		}
	});
}


// Options for tSNE plot //
function populateTSNEDropdowns() {
    var tsne_types_dropdown = $('#tsne_types');
    var tsne_cluster_dropdown = $('#tsne_clustering');

    $.getJSON('/tsne_options?q='+ensemble, function(data) {
        $.each(data['tsne_types'], function(key, val) {
            tsne_types_dropdown.append(
                $('<option></option>').val(val).text(val)
                );
        });
        $.each(data['cluster_types'], function(key, val) {
            tsne_cluster_dropdown.append(
                $('<option></option>').val(val).text(val)
                );
        });
    });
}


function randomizeClusterColors() {
    $('#randomizeColorBtn').click(function() {
        var grouping = $('#tsneGrouping option:selected').val();
        $.ajax("/plot/delete_cache/"+ensemble+"/"+grouping);
        if (grouping === 'biosample'){
            storeUpdate(empty=true);
            $.ajax({
                type: "GET",
                url: './plot/randomize_colors?n='+num_colors,
                success: function(data){
                    for(i = 0; i < data['num_colors']; i++){
                        var group = 'cluster_color_' + String(i);
                        var update = {
                            'marker.color': data['colors'][i]
                        };
                        Plotly.restyle("plot-2d-cluster", update, data[group]);
                        if($('#toggle-3d').prop('checked')) {
                            Plotly.restyle("plot-3d-cluster", update, data[group]);
                        }
                        else{
                            storeUpdate(update,data[group]);
                        }
                    }
                }
            });
        }
        else {
            loadClusterPlots();
        }
    });
}

function updateGeneElements(updateMCHCluster=true) {
    buttons = document.getElementsByClassName('modebar-btn');
    var geneSelected = $('#geneName option:selected').val();
    if (geneSelected != 'Select..' && $("#geneName").select2('data').length > 0) {
        $('#tSNE_cluster_div').addClass("col-md-6");
        $('#tSNE_cluster_div').removeClass("col-md-8 col-md-offset-2");
        $('#methyl_scatter_div, #methyl_graphs_div').show();
        $('#orthologsToggle').bootstrapToggle('off');
        $('#orthologsToggle').bootstrapToggle('enable');
        try{
            buttons[8].click();
        }
        catch(err) {
            console.log(err);
        }

        var lastViewedGenes = [];
        for(i=0; i<$('#geneName').select2('data').length; i++){
            lastViewedGenes.push({gene_name: $('#geneName option:selected')[i].text, gene_id: $('#geneName option:selected')[i].value});
        }
        if (typeof(Storage) !== 'undefined') {
            storage.save('lastViewedGenes', lastViewedGenes, 5);  // store last viewed genes for 5 minutes
        }
        if (updateMCHCluster){
            updateMCHClusterPlot();
        }
        if($("#geneName").select2('data').length > 1) {
            $('#normalize-heatmap').show();
            $('#normalize-toggle').prop('disabled', false);
            createHeatmap();
            updateDataTable("Select..");
            $('#epiBrowserLink').addClass('disabled');
        }
        else{
            $('#epiBrowserLink').removeClass('disabled');
            $('#normalize-heatmap').hide();
            $('#normalize-toggle').prop('disabled', true);
            updateOrthologToggle();
            updateMCHBoxPlot();
            updateDataTable($('#geneName option:selected').val());

            $.ajax({
                url: './gene/id/'+ensemble+'?q='+geneSelected,
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
    else{
        $('#methyl_scatter_div, #methyl_graphs_div').hide();
        $('#tSNE_cluster_div').addClass("col-md-8 col-md-offset-2");
        $('#tSNE_cluster_div').removeClass("col-md-6");
        try{buttons[8].click();}
        catch(e) {}
    }
}

function loadClusterPlots() {
    var grouping = $('#tsne_grouping').val();
    var clustering = $('#tsne_clustering').val();
    var tsne_type = $('#tsne_types').val();

    $.ajax({
        type: "GET",
        url: './plot/tsne/'+ensemble+'/'+tsne_type+'/'+grouping+'/'+clustering, 
        success: function(data) {
            num_colors = getMax(data["traces"], "legendgroup");
            Plotly.newPlot("plot-2d-cluster", Object.values(data["traces"]), data["layout"], {showLink: false});
            $('#loading_2dtsne').html("");
        }
    });
}

function updateMCHClusterPlot() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var pValues = pSlider.getValue();
    var genes = $("#geneName").select2('data');
    var clustering = $('#tsne_clustering').val();
    var tsne_type = $('#tsne_types').val();
    var genes_query = "";
    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    genes_query = genes_query.slice(0,-1);
    if ($('#geneName option:selected').val() != 'Select..') {
        $.ajax({
            type: "GET",
            url: './plot/scatter/'+ensemble+'/'+tsne_type+'/' +methylationType+ '/'+levelType+'/'+clustering+'/'+pValues[0]+'/'+pValues[1]+'?q='+genes_query,
            success: function(data) {
                $('#plot-mch-scatter').html("");
                $('#plot-mch-scatter').html(data);
            }
        });
    }
}

function updateOrthologToggle() {
    var geneSelected = $('#geneName option:selected').val();
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

function initDataTableClick() {
    $('#geneTable tbody').on('click', 'tr', function () {
        var id = $(this).attr('id');
        $.getJSON({
            url: './gene/id/'+ensemble+'?q='+id,
            success: function (data) {
                var option = new Option(data.gene_name, data.gene_id, true, true);
                var i;
                for(i=0; i < $("#geneName").select2('data').length; i++){
                    if($("#geneName").select2('data')[i].id === option.value){
                        return;
                    }
                }
                geneNameSelector.append(option);
                $('#epiBrowserLink').attr('href', generateBrowserURL(data));
                $('#epiBrowserLink').removeClass('disabled');
                updateGeneElements();
            }
        });
    });
}

function updateDataTable(geneSelected) {
    if (geneSelected !== 'Select..' || geneSelected !== "") {
        table = $('#geneTable').DataTable( {
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
                "dataSrc": ""
            },
            "rowId": 'gene_id',
            "columns": [
                { "data": "Rank" },
                { "data": "gene_name" },
                { "data": "Corr" },
            ],
        });
    }
    else {
        table.clear();
    }
}

function updateMCHBoxPlot() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var geneSelected = $('#geneName option:selected').val();
    var grouping = $('#tsne_grouping').val();
    var clustering = $('#tsne_clustering').val();
    if ($('#orthologsToggle').prop('checked')) {
        return updateMCHCombinedBoxPlot(mmu_gID, hsa_gID);
    }
    if ($('#outlierToggle').prop('checked')) {
        var outlierOption = 'outliers';
    } else {
        var outlierOption = 'false';
    }

    $.ajax({
        type: "GET",
        url: './plot/box/'+ensemble+'/'+methylationType+'/'+geneSelected+'/'+grouping+'/'+clustering+'/'+levelType+'/'+outlierOption,
        success: function(data) {
            $('#plot-mch-heat').html("");
            $('#mch_box_div').addClass("col-md-9");
            $('#gene_table_div').show();
            $('#plot-mch-box').html(data);
        }
    });

}

function updateMCHCombinedBoxPlot(mmu_gid, hsa_gid) {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    if ($('#outlierToggle').prop('checked')) {
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

function createHeatmap() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var pValues = pSlider.getValue();
    var genes = $("#geneName").select2('data');
    var genes_query = "";
    var grouping = $('#tsne_grouping').val();
    var clustering = $('#tsne_clustering').val();

    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    if ($('#normalize-toggle').prop('checked')) {
        var normalize = 'true';
    }
    else {
        var normalize = 'false';
    }
    genes_query = genes_query.slice(0,-1);

    $.ajax({
        type: "GET",
        url: './plot/heat/'+ensemble+'/'+methylationType+'/'+grouping+'/'+clustering+'/'+levelType+'/'+pValues[0]+'/'+pValues[1]+'?q='+genes_query+'&normalize='+normalize,
        success: function(data) {
            $('#plot-mch-box').html("");
            $('#gene_table_div').hide();
            $('#mch_box_div').removeClass("col-md-9");
            $('#plot-mch-heat').html(data);
            $('#outlierToggle').bootstrapToggle('disable');
        }
    });
}

function createHeatmapTwoSpecies() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var pValues = pSlider.getValue();
    var genes = $("#geneName").select2('data');
    var genes_query = "";
    
    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    if ($('#normalize-toggle').prop('checked')) {
        var normalize = 'true';
    }
    else {
        var normalize = 'false';
    }
    genes_query = genes_query.slice(0,-1);

    $.ajax({
        type: "GET",
        url: './plot/heat_two_ensemble/'+ensemble+'/'+methylationType+'/'+levelType+'/'+pValues[0]+'/'+pValues[1]+'?q='+genes_query+'&normalize='+normalize,
        success: function(data) {
            $('#plot-mch-box').html("");
            $('#gene_table_div').hide();
            $('#mch_box_div').removeClass("col-md-9");
            $('#plot-mch-heat').html(data);
            $('#outlierToggle').bootstrapToggle('disable');
        }
    });
}
