function initEnsembleDataTable() {

    var table = $('#ensemble_table').DataTable({
        "ajax": {
            "url": "/content/ensemble_list",
            "dataSrc": ""
        },
        "columns": [
            {
                "className": 'details-control dt-center',
                "orderable": false,
                "data": null,
                "defaultContent": '<i class="glyphicon glyphicon-plus-sign"></i>'
            },
            {"data": "ensemble_id"},
            {"data": "ensemble_name"},
            {"data": "total_methylation_cells"},
            {"data": "total_snATAC_cells"},
            {"data": "ABA_regions_acronym"},
            {"data": "slices"},
            {"data": "num_datasets"},
            {
                "className": 'redirect-control dt-center',
                "orderable": false,
                "data": "ensemble_name",
                "fnCreatedCell": function (nTd, sData, oData, iRow, iCol) {
                    $(nTd).html('<a href="/'+oData.ensemble_name+'"><i class="glyphicon glyphicon-eye-open"></i></a>');
                }
            },
        ]
    });

    $('#ensemble_table tbody').on('click', 'td.details-control', function() {
        var tr = $(this).closest('tr');
        var row = table.row(tr);

        if (row.child.isShown()) {
            row.child.hide();
            tr.removeClass('shown');
            $(this).html('<i class="glyphicon glyphicon-plus-sign"></i>');
        }
        else {
            row.child(format(row.data())).show();
            tr.addClass('shown');
            $(this).html('<i class="glyphicon glyphicon-minus-sign"></i>');
        }
    });
}

function format ( d ) {

    var snATAC_datasets = d.snATAC_datasets;
    if (d.snATAC_datasets === "") {
        snATAC_datasets = "None";
    }

    return (
        '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">'+
            '<tr>'+
                '<td>Datasets in ensemble (snmC-seq):</td>'+
                '<td>'+d.datasets+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td>Datasets in ensemble (snATAC-seq):</td>'+
                '<td>'+snATAC_datasets+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td>Brain regions:</td>'+
                '<td>'+d.ABA_regions_description+'</td>'+
            '</tr>'+
        '</table>'
    )
    
    
}
