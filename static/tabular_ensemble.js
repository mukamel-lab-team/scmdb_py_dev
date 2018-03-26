function initEnsembleDataTable() {

    var table = $('#ensemble-table').DataTable({
        "order": [[3, 'desc']], //Initially sort by "Total Cells (snmC-seq)" in descending order. 
        "pageLength": 25,
        "ajax": {
            "url": "/content/ensemble_list",
            "dataSrc": ""
        },
        "columns": [
            {
                "data": null,
                "className": 'details-control dt-center',
                "orderable": false,
                "defaultContent": '<i class="glyphicon glyphicon-plus-sign"></i>'
            },
            {"data": "ensemble_id"},
            {"data": "ensemble_name"},
            {"data": "total_methylation_cells", "className": 'dt-center'},
            {"data": "total_snATAC_cells", "className": 'dt-center'},
            {"data": "ABA_regions_acronym", "className": 'dt-center'},
            {"data": "slices", "className": 'dt-center'},
            {"data": "num_datasets", "className": 'dt-center'},
            {
                "data": "public_access_icon",
                "className": 'dt-center',
                "fnCreatedCell": function(nTd, sData, oData, iRow, iCol) {
                    $(nTd).html('<i class="'+oData.public_access_icon+'" style="color:'+oData.public_access_color+';"></i>');
                }
            },
            {
                "data": "ensemble_name",
                "className": 'redirect-control dt-center',
                "orderable": false,
                "fnCreatedCell": function (nTd, sData, oData, iRow, iCol) {
                    $(nTd).html('<a href="/'+oData.ensemble_name+'"><i class="glyphicon glyphicon-eye-open"></i></a>');
                }
            },
        ]
    });

    $('#ensemble-table tbody').on('click', 'td.details-control', function() {
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
                '<td><b>Datasets in ensemble (snmC-seq):</b></td>'+
                '<td>'+d.datasets+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Datasets in ensemble (snATAC-seq):</b></td>'+
                '<td>'+snATAC_datasets+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Brain region(s):</b></td>'+
                '<td>'+d.ABA_regions_description+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>View Dissection Drawings:</b></td>'+
                '<td><a href="https://drive.google.com/file/d/1dAUzB1GtKMUgmf_AInAGgI6OlgefUHok/preview" target="_blank">Link</a> (go to page '+d.slices+')</td>'+
            '</tr>'+
        '</table>'
    )
}
