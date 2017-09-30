function loadClusterPlot() {
    var grouping = $('#tsneGrouping option:selected').val();
    $.ajax({
        type: "GET",
        url: './plot/cluster/' + species + '/' + grouping,
        success: function(data) {
            $('#plot-cluster').html(data);
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

    // Initialise selector
    if (typeof(Storage) !== 'undefined') {
        var defaultGene = localStorage.getItem('lastViewed');
        if (defaultGene == null) {
            defaultGene = 'Gad2'; // No entry, set default to Gad2.
        }
    } else {
        // Browser has no localStorage support, we'll just do Gad2.
        var defaultGene = 'Gad2';
    }

    $.getJSON({
        url: './gene/names/' + species + '?q=' + defaultGene,
        success: function(data) {
            data.forEach(function(gene) {
                var option = new Option(gene.geneName, gene.geneID, true, true);
                geneNameSelector.append(option);
                $('#epiBrowserLink').attr('href', generateBrowserURL(gene));
                $('#epiBrowserLink').removeClass('disabled');
            });
            geneNameSelector.trigger('change');
        }
    });
}

function updateMCHClusterPlot() {
    var levelType = $('input[name=levels]').filter(':checked').val();
    var geneSelected = $('#geneName option:selected').val();
    var pValues = pSlider.getValue();
    if (geneSelected != 'Select..') {
        $.ajax({
            type: "GET",
            url: './plot/mch/' + species + '/' + geneSelected + '/' + levelType + '/' + pValues[0] + '/' + pValues[1],
            success: function(data) {
                $('#plot-mch-scatter').html(data);
            }
        });
    }
}


function updateMCHBoxPlot() {
    var levelType = $('input[name=levels]').filter(':checked').val();
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
        url: './plot/box/' + species + '/' + geneSelected + '/' + levelType + '/' + outlierOption,
        success: function(data) {
            $('#plot-mch-box').html(data);
        }
    });

}

function updateMCHCombinedBoxPlot(mmu_gid, hsa_gid) {
    var levelType = $('input[name=levels]').filter(':checked').val();
    if ($('#outlierToggle').prop('checked')) {
        var outlierOption = 'outliers';
    } else {
        var outlierOption = 'false';
    }
    $.ajax({
        type: "GET",
        url: './plot/box_combined/' + mmu_gid + '/' + hsa_gid + '/' + levelType + '/' + outlierOption,
        success: function(data) {
            $('#plot-mch-box').html(data);
        }
    });

}

function updateGeneElements() {
    var geneSelected = $('#geneName option:selected').val();
    if (geneSelected != 'Select..') {
        if (typeof(Storage) !== 'undefined') {
            localStorage.setItem('lastViewed', $('#geneName option:selected').text());
        }
        updateMCHClusterPlot();
        updateOrthologToggle();
        updateMCHBoxPlot();
        try {
            var annojURL;
            geneSearchCache.forEach(function(gene) {
                if (gene.geneID == geneSelected) {
                    annojURL = generateBrowserURL(gene);
                }
            });
            $('#epiBrowserLink').attr('href', annojURL);
        } catch (e) {}

    }
}

function generateBrowserURL(gene) {
    if (species == 'mmu') {
        var base = 'http://brainome.ucsd.edu/annoj/CEMBA/index_mm.html'; // Mouse
    } else {
        var base = 'http://brainome.ucsd.edu/annoj/sc_wgbs/index_hs.html'; // Human
    }

    if (gene.strand == '+') {
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
            if (data.mmu_gID == null || data.hsa_gID == null) {
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
            url: './gene/id/' + species + '?q=' + id,
            success: function (data) {
                var option = new Option(data.geneName, data.geneID, true, true);
                geneNameSelector.append(option);
                $('#epiBrowserLink').attr('href', generateBrowserURL(data));
                $('#epiBrowserLink').removeClass('disabled')
                updateGeneElements();
                updateDataTable();
            }
        });
    });
}

function updateDataTable() {
    var geneSelected = $('#geneName option:selected').val();
    if (geneSelected != 'Select..') {
        var table = $('#geneTable').DataTable( {
            "destroy": true,
            "processing": true,
            "ordering": false,
            "lengthChange": false,
            "dom": "<'row'<'col-md-3 'f>>" +
                    "<'row'<'col-md-3'tr>>" +
                    "<'row'<'col-md-3'i>>" + 
                    "<'row'<'col-md-3'p>>",
            "pagingType": "simple",
            "ajax": {
                "url": "./gene/corr/" + species + "/" + geneSelected,
                "dataSrc": ""
            },
            rowId: 'geneID',
            "columns": [
                { "data": "Rank" },
                { "data": "geneName" },
                { "data": "Corr" },
            ]
        });
    }
}