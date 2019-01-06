function initEnsembleDataTable() {
    
    const numberWithCommas = (x) => {
      return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
    }
    
    var table = $('#ensemble-table').DataTable({
        "order": [[3, 'desc']], //Initially sort by "Total Cells (snmC-seq)" in descending order. 
        "pageLength": 100,
        "ajax": {
            //"url": "/portal/content/ensembles",
            "url": $SCRIPT_ROOT+'/content/ensembles?region=' + region,
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
            {
                "data": "ensemble_name",
                "className": 'redirect-control dt-center',
                "orderable": false,
                "fnCreatedCell": function (nTd, sData, oData, iRow, iCol) {
                    if ( oData.annoj_exists == 1 )   {
                        $(nTd).html('<a href="https://brainome.ucsd.edu/annoj_private/CEMBA/cemba.php?ens='+oData.ensemble_id+'" target="_blank"><i class="fas fa-eye"></i></a>');
                    } else {
                        $(nTd).html('<i class="fas fa-eye-slash" style="color: gray"></i>')
                    }
                }
            },
        ],
        "fixedHeader": {
            "header": false,
            "footer": false,
        },
        "footerCallback": function ( row, data, start, end, display ) {
            var api = this.api(), data;

            // Total mC cells
            grand_total_methylation_cells = api
                .column( 3 )
                .data()
                .reduce( function (a, b) {
                    return a+b;
                }, 0 );
 
            // Total snATAC cells
            grand_total_snATAC_cells = api
                .column( 4 )
                .data()
                .reduce( function (a, b) {
                    return a+b;
                }, 0 );
 
            grand_total_ensembles = api
                .column( 1 )
                .data().map(function(e){ return (e>0); })
                .reduce( function (a, b) {
                    return a+b;
                }, 0 );
 
            // Update footer
            $( api.column( 10 ).footer() ).html(
                numberWithCommas(grand_total_methylation_cells) +' mC cells, '+ 
                numberWithCommas(grand_total_snATAC_cells) +
                ' ATAC cells across '+grand_total_ensembles+' ensembles'
            );
        }
    } );
    

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

    var out = '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">'+
            '<tr>'+
                '<td><b>Description:</b></td>'+
                '<td>'+d.description+'</td>'+
            '</tr>'+
            '<tr>'+
                '<td><b>Datasets in ensemble (snmC-seq):</b></td>'+
                '<td>'+
                    '<table>';
    if (d.datasets_rs1 != "") { out = out+'<tr><td><b>RS1</b></td><td>'+d.datasets_rs1+'</td></tr>' }
    if (d.datasets_rs2 != "") { out = out+'<tr><td><b>RS2</b></td><td>'+d.datasets_rs2+'</td></tr>' }
    out = out + 
                    '</table>'+
                '</td>'+
            '</tr>'
    if (d.snATAC_datasets_rs1+d.snATAC_datasets_rs2 != "") { out = out+
            '<tr>'+
                '<td><b>Datasets in ensemble (snATAC-seq):</b></td>'+
                '<td>'+
                    '<table>'}
    if (d.snATAC_datasets_rs1 != "") { out = out+'<tr><td><b>RS1</b></td><td>'+d.snATAC_datasets_rs1+'</td></tr>' }
    if (d.snATAC_datasets_rs2 != "") { out = out+'<tr><td><b>RS2</b></td><td>'+d.snATAC_datasets_rs2+'</td></tr>' }
    if (d.snATAC_datasets_rs1+d.snATAC_datasets_rs2 != "") {
            out = out + 
                        '</table>'+
                '</td>'+
            '</tr>'}
            out = out + 
            '<tr>'+
                '<td><b>Dissection region(s):</b></td>'+
                '<td>'+d.ABA_regions_description+'</td>'+
            '</tr>'
    if (d.target_regions_rs2_descriptive != "") {
        out = out+
            '<tr>'+
                '<td><b>Target region(s):</b></td>'+
                '<td>'+d.target_regions_rs2_descriptive+'</td>'+
            '</tr>'
        }
        out = out + 
            '<tr>'+
                '<td><b>View Dissection Drawings:</b></td>'+
                '<td><a href="https://drive.google.com/file/d/1dAUzB1GtKMUgmf_AInAGgI6OlgefUHok/preview" target="_blank">Link</a> (go to page '+d.slices+')</td>'+
            '</tr>'+
        '</table>'
    return out
}
