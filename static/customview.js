var trace_3d, layout_3d;
var updates_3d = [];
var groups_3d = [];

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
	}
}

function getInitGeneInfo(species, geneID) {
    return $.getJSON('./gene/id/' + species + '?q=' + geneID).then(function(data) {
        return {
            geneName: data.geneName,
            geneID: data.geneID
        }
    });
}

function loadClusterPlots() {
    var grouping = $('#tsneGrouping option:selected').val();
    $.ajax({
        type: "GET",
        url: './plot/cluster/' + species + '/' + grouping,
        success: function(data) {
            Plotly.newPlot("plot-2d-cluster", Object.values(data["traces_2d"]), data["layout2d"]);
            $('#loading_2dtsne').html("");
            if("traces_3d" in data){
                if($('#toggle-3d').prop('checked')) {
                    Plotly.newPlot("plot-3d-cluster", Object.values(data["traces_3d"]), data["layout3d"]);
                    save3DData(data["traces_3d"],data["layout3d"]);
                }
                else {
                    save3DData(data["traces_3d"], data["layout3d"]);
                    $('#toggle-3d').removeAttr("disabled");
                }
            }
            else {
                save3DData(null, null);
                $('#toggle-3d').attr("disabled", true);
            }
        }
    });
}

function initGeneNameSearch() {
    geneNameSelector = $('#geneName').select2({
        placeholder: 'Search..',
        allowClear: true,
        ajax: {
            url: './gene/names/' + species,
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
                            text: gene.geneName,
                            id: gene.geneID
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
    if (!defaultGene) {
        //no entry or browser does not support localStorage, set default to none
        defaultGene = [];
    }
    if(defaultGene !== []){
        var numGenes = defaultGene.length;
        for (i = 0; i < numGenes; i++) {
            $.ajax({
                url: './gene/id/' + species + '?q=' + defaultGene[i].geneID,
                dataType: 'json',
                async: false,
                success: function(data) {
                    console.log(data);
                    if (typeof(data.geneName) !== 'undefined' && typeof(data.geneID) !== 'undefined') {
                        var option = new Option(data.geneName, data.geneID, true, true);
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
		url: './gene/modules/' + species,
		success: function(data){
			data.forEach(function(gene) {
				var option = new Option(gene.module, gene.module, false, false)
				geneModuleSelector.append(option);
			});
		}
	});
}

function updateSearchWithModules(module) {
	$.getJSON({
		url: './gene/modules/' + species + '?q=' + module.id,
		success: function (data) {
			data.forEach(function(gene) {
				var option = new Option(gene.geneName, gene.geneID, true, true);
				geneNameSelector.append(option);
			});
		}
	});
}

function updateMCHClusterPlot() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var pValues = pSlider.getValue();
    var genes = $("#geneName").select2('data');
    var genes_query = "";
    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    genes_query = genes_query.slice(0,-1);
    if ($('#geneName option:selected').val() != 'Select..') {
        $.ajax({
            type: "GET",
            url: './plot/scatter/' + species + '/' + methylationType +  '/' + levelType + '/' + pValues[0] + '/' + pValues[1] + '?q=' + genes_query,
            success: function(data) {
                $('#plot-mch-scatter').html("");
                $('#plot-mch-scatter').html(data);
            }
        });
    }
}


function updateMCHBoxPlot() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var geneSelected = $('#geneName option:selected').val();
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
        url: './plot/box/' + species + '/' + methylationType + '/' + geneSelected + '/' + levelType + '/' + outlierOption,
        success: function(data) {
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
        url: './plot/box_combined/' + species + '/' + methylationType + '/' + mmu_gid + '/' + hsa_gid + '/' + levelType + '/' + outlierOption,
        success: function(data) {
            $('#plot-mch-box').html(data);
        }
    });
}

function updateGeneElements() {
    buttons = document.getElementsByClassName('modebar-btn');
    var geneSelected = $('#geneName option:selected').val();
    if (geneSelected != 'Select..' && $("#geneName").select2('data').length > 0) {
        $('#tSNE_cluster_div').addClass("col-md-6");
        $('#tSNE_cluster_div').removeClass("col-md-8 col-md-offset-2");
        $('#methyl_scatter_div, #methyl_graphs_div').show();
        try{
            buttons[8].click();
        }
        catch(err) {
            console.log(err);
        }

        var lastViewedGenes = [];
        for(i=0; i<$('#geneName').select2('data').length; i++){
            lastViewedGenes.push({geneName: $('#geneName option:selected')[i].text, geneID: $('#geneName option:selected')[i].value});
        }
        if (typeof(Storage) !== 'undefined') {
            storage.save('lastViewedGenes', lastViewedGenes, 5);  // store last viewed genes for 5 minutes
        }
        updateMCHClusterPlot();
        console.log($("#geneName").select2('data').length); 
        if($("#geneName").select2('data').length > 1) {
            createHeatMap();
            updateDataTable("Select..");
            $('#epiBrowserLink').addClass('disabled');
        }
        else{
            $('#outlierToggle').bootstrapToggle('enable');
            $('#orthologsToggle').bootstrapToggle('enable');
            $('#epiBrowserLink').removeClass('disabled');
            updateOrthologToggle();
            updateMCHBoxPlot();
            updateDataTable($('#geneName option:selected').val());
        }
        try {
            var annojURL;
            geneSearchCache.forEach(function(gene) {
                if (gene.geneID === geneSelected) {
                    annojURL = generateBrowserURL(gene);
                }
            });
            $('#epiBrowserLink').attr('href', annojURL);
        } catch (e) {}
    }
    else{
        $('#methyl_scatter_div, #methyl_graphs_div').hide();
        $('#tSNE_cluster_div').addClass("col-md-8 col-md-offset-2");
        $('#tSNE_cluster_div').removeClass("col-md-6");
        try{buttons[8].click();}
        catch(e) {}
    }
}

function generateBrowserURL(gene) {
    if (species === 'mmu') {
        var base = 'http://brainome.ucsd.edu/annoj/CEMBA/index_mm.html'; // Mouse
    } else {
        var base = 'http://brainome.ucsd.edu/annoj/sc_wgbs/index_hs.html'; // Human
    }

    if (gene.strand === '+') {
        var position = gene.start;
    } else {
        var position = gene.end;
    }
    return base + '?assembly=' + gene.chrom + '&position=' + position;
}

function updateOrthologToggle() {
    var geneSelected = $('#geneName option:selected').val();
    $.ajax({
        type: "GET",
        url: './gene/orthologs/' + species + '/' + geneSelected,
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

function save3DData(trace, layout){
    trace_3d = trace;
    layout_3d = layout;
}

function display3DPlotToggle() {
    if ($('#toggle-3d').prop('checked')){
        $('#loading-3d-plot').html("loading..");
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

function initDataTableClick() {
    $('#geneTable tbody').on('click', 'tr', function () {
        var id = $(this).attr('id');
        $.getJSON({
            url: './gene/id/' + species + '?q=' + id,
            success: function (data) {
                var option = new Option(data.geneName, data.geneID, true, true);
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
                "url": "./gene/corr/" + species + "/" + geneSelected,
                "dataSrc": ""
            },
            "rowId": 'geneID',
            "columns": [
                { "data": "Rank" },
                { "data": "geneName" },
                { "data": "Corr" },
            ],
        });
    }
    else {
        table.clear();
    }
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

function randomizeClusterColors() {
    $('#randomizeColorBtn').click(function() {
        var grouping = $('#tsneGrouping option:selected').val();
        if (grouping === 'biosample'){
            storeUpdate(empty=true);
            $.ajax({
                type: "GET",
                url: './plot/randomize_colors',
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
            $.get("/plot/delete_cache/" + species + "/" + grouping);
            loadClusterPlots();
        }
    });
}

function createHeatMap() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var methylationType = $('input[name=mType]').filter(':checked').val();
    var pValues = pSlider.getValue();
    var genes = $("#geneName").select2('data');
    var genes_query = "";
    for (i = 0; i < genes.length; i++) {
        genes_query += (genes[i].id + "+");
    }
    genes_query = genes_query.slice(0,-1);
    $.ajax({
        type: "GET",
        url: './plot/heat/' + species + '/' + methylationType + '/' + levelType + '/' + pValues[0] + '/' + pValues[1] + '?q=' + genes_query,
        success: function(data) {
            $('#plot-mch-box').html(data);
            $('#outlierToggle').bootstrapToggle('disable');
            $('#orthologsToggle').bootstrapToggle('disable');
        }
    });
}
