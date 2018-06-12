function initEnsembleDataTable() {

    var table = $('#ensemble-table').DataTable({
        "order": [[3, 'desc']], //Initially sort by "Total Cells (snmC-seq)" in descending order. 
        "pageLength": 25,
        "ajax": {
            //"url": "/portal/content/ensembles",
            "url": $SCRIPT_ROOT+'/content/ensembles',
            "dataSrc": ""
        },
        "columns": [
            {
                "data": null,
                "className": 'details-control dt-center',
                "orderable": false,
                "defaultContent": '<i class="fas fa-plus-circle"></i>'
            },
            {"data": "ensemble_id"},
            {"data": "ensemble_name"},
            {"data": "total_methylation_cells", "className": 'dt-center'},
            {"data": "total_snATAC_cells", "className": 'dt-center'},
            {"data": "ABA_regions_acronym", "className": 'dt-center'},
            {"data": "slices", "className": 'dt-center'},
            {"data": "num_datasets", "className": 'dt-center'},
            {"data": "target_regions_rs2_acronym", "className": 'dt-center'},
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
                    $(nTd).html('<a href="'+$SCRIPT_ROOT+'/'+oData.ensemble_name+'"><i class="fas fa-eye"></i></a>');
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
            $(this).html('<i class="fas fa-plus-circle"></i>');
        }
        else {
            row.child(format(row.data())).show();
            tr.addClass('shown');
            $(this).html('<i class="fas fa-minus-circle"></i>');
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
                '<td><b>Description:</b></td>'+
                '<td>'+d.description+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Datasets in ensemble (snmC-seq):</b></td>'+
                '<td>'+
                    '<table>'+
                        '<tr><td><b>RS1</b></td><td>'+d.datasets_rs1+'</td></tr>'+
                        '<tr><td><b>RS2</b></td><td>'+d.datasets_rs2+'</td></tr>'+
                    '</table>'+
                '</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Datasets in ensemble (snATAC-seq):</b></td>'+
                '<td>'+
                    '<table>'+
                        '<tr><td><b>RS1</b></td><td>'+d.snATAC_datasets_rs1+'</td></tr>'+
                        '<tr><td><b>RS2</b></td><td>'+d.snATAC_datasets_rs2+'</td></tr>'+
                    '</table>'+
                '</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Dissection region(s):</b></td>'+
                '<td>'+d.ABA_regions_description+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Target region(s):</b></td>'+
                '<td>'+d.target_regions_rs2_descriptive+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>View Dissection Drawings:</b></td>'+
                '<td><a href="https://drive.google.com/file/d/1dAUzB1GtKMUgmf_AInAGgI6OlgefUHok/preview" target="_blank">Link</a> (go to page '+d.slices+')</td>'+
            '</tr>'+
        '</table>'
    )
}
