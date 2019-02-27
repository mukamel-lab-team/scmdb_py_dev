function initEnsembleDataTable() {
    
    const numberWithCommas = (x) => {
      return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
    }
    
    var table = $('#ensemble-table').DataTable({
        "order": [[3, 'desc']], //Initially sort by "Total Cells (snmC-seq)" in descending order. 
        "pageLength": 50,
        "ajax": {
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
                "data": "annoj_exists",
                "className": 'dt-center',
                "fnCreatedCell": function (nTd, sData, oData, iRow, iCol) {
                    if ( oData.annoj_exists == 1 )  {
                        $(nTd).html('<a href="https://brainome.ucsd.edu/annoj_private/CEMBA/cemba.php?ens='+oData.ensemble_id+'" target="_blank"><i class="fas fa-dna"></i></a>');
                    } else {
                        $(nTd).html('<i class="fas fa-dna" style="color: gray"></i>');
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
                grand_total_ensembles+' ensembles including '+numberWithCommas(grand_total_methylation_cells) +' mC cells, '+ 
                numberWithCommas(grand_total_snATAC_cells) +
                ' ATAC cells across (NOTE some cells appear in multiple ensembles)'
            );
        }
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
    var slices=[... new Set(d.slices.split(', ').map(x => x.charAt(0)))] // Get the unique slices contributing to this ensemble
    out = out + 
        '<tr>'+
            '<td><a target=”_blank" href="https://brainome.ucsd.edu/CEMBA_dissection_images/index.php?slideIndex='+Math.min.apply(Math, slices)+'">'+
            'Dissection drawings for all slices</a></td>'+
        '</tr>'+
        '<tr>'
    var url=''
    for (i=0; i<slices.length; i++) {
        url='https://brainome.ucsd.edu/CEMBA_dissection_images/CEMBA_Slice'+slices[i]
        out = out+'<td> '+
          '<a href="https://brainome.ucsd.edu/CEMBA_dissection_images/index.php?slideIndex='+slices[i]+'" target="_blank">'+
          '<img src="'+url+'_sm.png" ' +
          ' width=200 alt="CEMBA Slice'+slices[i]+'"></a> </td>'
    }
    out = out + '</tr></table>'
    return out
}
