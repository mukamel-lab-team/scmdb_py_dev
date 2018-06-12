var table; 

function initDataTable() {
    
    table = $('#dataset-table').DataTable( {
        "stateSave": true,
        "select": true,
        "order": [[2, 'Desc']], //Initially sort by date added.
        "pageLength": 25,
        "ajax": {
            "url": "/content/datasets/all",
            "dataSrc": ""
        },
        "columns": [
            {
                "data": null,
                "className": 'details-control dt-center',
                "orderable": false,
                "defaultContent": '<i class="fas fa-plus-circle"></i>'
            },
            {"data": "dataset_name"},
            {"data": "research_segment", "className": 'dt-center'},
            {"data": "sex", "className": 'dt-center'},
            {"data": "methylation_cell_count", "className": 'dt-center'},
            {"data": "snATAC_cell_count", "className": 'dt-center'},
            {"data": "ABA_regions_acronym", "className": 'dt-center'},
            {"data": "slice", "className": 'dt-center'},
            {"data": "target_region_acronym", "className": 'dt-center'},
            {"data": "date_added", "className": 'dt-center'},
        ]
    } );

    $('#dataset-table tbody').on('click', 'td.details-control', function() {
        var tr = $(this).closest('tr');
        var row = table.row(tr);

        if (row.child.isShown()) {
            row.child.hide();
            tr.removeClass('shown');
            $(this).html('<i class="fas fa-plus-circle"></i>');
        }
        else {
            row.child(formatChildRows(row.data())).show();
            tr.addClass('shown');
            $(this).html('<i class="fas fa-minus-circle"></i>');
        }
    });
    return table;
}

function formatChildRows ( d ) {

    var snATAC_datasets = d.snATAC_datasets;
    if (d.snATAC_datasets === "") {
        snATAC_datasets = "None";
    }

    return (
        '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">'+
            '<tr>'+
                '<td><b>Description:</b></td>'+
                '<td>'+d.description+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Dissection region(s):</b></td>'+
                '<td>'+d.ABA_regions_descriptive+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Virus target region(s):</b></td>'+
                '<td>'+d.target_region_descriptive+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>View Dissection Drawings:</b></td>'+
                '<td><a href="https://drive.google.com/file/d/1dAUzB1GtKMUgmf_AInAGgI6OlgefUHok/preview" target="_blank">Link</a> (go to page '+d.slice+')</td>'+
            '</tr>'+
        '</table>'
    )
}
        
function enableDataTableSelect() {

    table.on('select', function( e, dt, type, indexes) {
        let selected_datasets = [];
        let selected_slices_set = new Set();
        let selected_brain_regions_set = new Set();
        let selected_target_regions_set = new Set();
        let total_methylation_cells = 0;
        let total_snATAC_cells = 0;

        let data = table.rows({selected:true}).data();
        for (i = 0; i < data.length; i++) {
            selected_datasets.push(data[i]['dataset_name']);
            selected_slices_set.add(data[i]['slice']);
            selected_brain_regions_set.add(data[i]['ABA_regions_acronym']);
            selected_target_regions_set.add(data[i]['target_region_acronym']);
            total_methylation_cells += data[i]['methylation_cell_count'];
            total_snATAC_cells += data[i]['snATAC_cell_count'];
        }
        $("#num-datasets").text("("+data.length+" selected)");
        $("#selected-datasets").val(selected_datasets.join(", "));
        $("#selected-datasets-slices").val([...selected_slices_set].join(", "));
        $("#selected-datasets-brain-regions").val([...selected_brain_regions_set].join(", "));
        $("#selected-datasets-target-regions").val([...selected_target_regions_set].join(", "));
        $("#selected-datasets-total-methylation-cells").val(total_methylation_cells);
        $("#selected-datasets-total-snATAC-cells").val(total_snATAC_cells);

        $("#check-ensemble").show();
        $("#check-button").prop("disabled", false);
        $("#check-button-text").show();
        $("#verification-message").hide();
        $("#final-submit-button").prop("disabled", true);
        $("#final-submit-button").hide();
    });

    table.on('deselect', function( e, dt, type, indexes) {
        let selected_datasets = [];
        let selected_slices_set = new Set();
        let selected_brain_regions_set = new Set();
        let selected_target_regions_set = new Set();
        let total_methylation_cells = 0;
        let total_snATAC_cells = 0;

        let data = table.rows({selected:true}).data();
        for (i = 0; i < data.length; i++) {
            selected_datasets.push(data[i]['dataset_name']);
            selected_slices_set.add(data[i]['slice']);
            selected_brain_regions_set.add(data[i]['ABA_regions_acronym']);
            selected_target_regions_set.add(data[i]['target_region_acronym']);
            total_methylation_cells += data[i]['methylation_cell_count'];
            total_snATAC_cells += data[i]['snATAC_cell_count'];
        }
        $("#num-datasets").text("("+data.length+" selected)");
        $("#selected-datasets").val(selected_datasets.join(", "));
        $("#selected-datasets-slices").val([...selected_slices_set].join(", "));
        $("#selected-datasets-brain-regions").val([...selected_brain_regions_set].join(", "));
        $("#selected-datasets-target-regions").val([...selected_target_regions_set].join(", "));
        $("#selected-datasets-total-methylation-cells").val(total_methylation_cells);
        $("#selected-datasets-total-snATAC-cells").val(total_snATAC_cells);
    });
}
        

function verifySimilarity() {
    $("#check-button-text").hide();
    $("#check-button-loader").show();
    $("#please-wait-message").show();

    $.getJSON({
        url: '/content/check_ensembles/'+$("#ensemble-name-input").val().trim()+'/'+$("#selected-datasets").val().replace(/, /g, '+'),
        success: function(data) {
            $("#please-wait-message").hide();
            $("#check-button-loader").hide();

            if (data['result'] === 'success') {
                $("#check-ensemble").hide();
                $("#verification-message").removeClass("warning error success").addClass("info");
                $("#verification-header").text("New ensemble is ready for creation!");
                $("#final-submit-button").prop("disabled", false);
                $("#final-submit-button").show()
            }
            else if (data['result'] === 'warning'){
                $("#verification-message").removeClass("error info success").addClass("warning");
                $("#verification-header").text("Warning");
                $("#final-submit-button").prop("disabled", false);
                $("#final-submit-button").show()
            }
            else {
                $("#check-button-text").show();
                $("#check-button").prop("disabled", false);
                $("#verification-message").removeClass("warning info success").addClass("error");
                $("#verification-header").text("Please try again");
            }

            $("#verification-message-text").text(data['reason']);
            $("#verification-message").show();
        }
    });
}

$(".ui.form").form( {
    on: 'blur',
    onSuccess: function() {
        $("#check-button").prop("disabled", true);
        verifySimilarity();
    },
    fields: {
        ensembleName: {
            identifier: 'ensemblename',
            rules: [{
                type: 'regExp[/^CEMBA_\\w+$/]',
                prompt: 'Invalid ensemble name format.',
                //value: /^CEMBA_\w*/,
            }]
        },
        ensembleDescription: {
            identifier: 'ensemble-description',
            rules: [{
                type: 'maxLength[255]',
                prompt: 'Description cannot exceed 255 characters.',
            }]
        },
        selectedDatasets: {
            identifier: 'selected-datasets',
            rules: [{
                type: 'empty',
                prompt: 'Must selected at least one dataset.',
            }]
        },
    }
});

$("#final-submit-button").on('click', function() {
    console.log($("#ensemble-description").val());
    $.getJSON({
        url: '/submit_new_ensemble/'+$("#ensemble-name-input").val().trim()+'/'+$("#selected-datasets").val().replace(/, /g, '+')+'?description='+$("#ensemble-description").val(),
        success: function(data) {
            if (data['result'] === true) {
                $("#verification-message").removeClass("warning info error").addClass("success");
                $("#verification-header").text("Success!");
                $("#verification-message-text").html("Your request has been sent. The data for <b>"+$("#ensemble-name-input").val()+"</b> will be available within the next 48 hours. ");
                $("#verification-message").show();
                $("#final-submit-button").prop('disabled', true);
            }
            else {
                $("#verification-message").removeClass("warning info success").addClass("error");
                $("#verification-header").text("Something went wrong!");
                $("#verification-message-text").text("A server-side error has occured and we were unable to submit your request. Please try again later.");
                $("#verification-message").show();
            }
        }
    });
});

$("#ensemble-description").on("change keyup paste", function() {
    var currentVal = $(this).val();
    if (currentVal == oldVal) {
        return; //check to prevent multiple simultaneous triggers
    }
    var oldVal = currentVal;

    const maxChars = 255;
    let charsRemaining = maxChars - $(this).val().length;
    $("#ensemble-description-chars-remaining").text(charsRemaining);

});

$("#ensemble-name-input").on("change keyup paste", function() {
    var currentVal = $(this).val();
    if (currentVal == oldVal) {
        return; //check to prevent multiple simultaneous triggers
    }
    var oldVal = currentVal;

    $("#check-ensemble").show();
    $("#check-button").prop("disabled", false);
    $("#check-button-text").show();
    $("#verification-message").hide();
    $("#final-submit-button").prop("disabled", true);
    $("#final-submit-button").hide();
});
