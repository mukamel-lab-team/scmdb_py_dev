function initDatasetDataTableRS1() {
    
    const numberWithCommas = (x) => {
      return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
    }

    var table = $('#dataset-table-rs1').DataTable( {
        "order": [[7, 'desc']], //Initially sort by date added.
        "pageLength": 50,
        "ajax": {
            "url": $SCRIPT_ROOT+"/content/datasets/rs1",
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
            {"data": "sex", "className": 'dt-center'},
            {"data": "methylation_cell_count", "className": 'dt-center'},
            {"data": "snATAC_cell_count", "className": 'dt-center'},
            {"data": "ABA_regions_acronym", "className": 'dt-center'},
            {"data": "slice", "className": 'dt-center'},
            {"data": "date_added", "className": 'dt-center'},
        ],      
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
 
            // Update footer
            $( api.column( 7 ).footer() ).html(
                numberWithCommas(grand_total_methylation_cells) +' mC cells, '+ 
                numberWithCommas(grand_total_snATAC_cells) +
                ' ATAC cells'
            );
        }
    });

    $('#dataset-table-rs1 tbody').on('click', 'td.details-control', function() {
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
    
    var slices=[... new Set(d.slice.split(', ').map(x => x.charAt(0)))] // Get the unique slices contributing to this ensemble
    out = '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">'+
        '<tr>'+
            '<td><b>Description:</b></td>'+
            '<td>'+d.description+'</td>'+
        '</tr>'+
        '<tr>'+
            '<td><b>Dissection region(s):</b></td>'+
            '<td>'+d.ABA_regions_descriptive+'</td>'+
        '</tr>'+
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
        
