/* rowGetter, label, dataKey, headerRenderer */
import React from 'react';
import {Table, Column, Cell} from 'fixed-data-table-2';
import "fixed-data-table-2/dist/fixed-data-table.css";
import ReactTooltip from 'react-tooltip';
import {findDOMNode} from 'react-dom';


class MyTextCell extends React.Component {
    render() {
        const {rowIndex, field, data} = this.props;
        return (
            <Cell {...this.props}>
                {data[rowIndex][field]}
            </Cell>
       );
    }
}

class MyLinkCell extends React.Component {
    render() {
        const {rowIndex, field, data} = this.props;
        const link = '/' + data[rowIndex][field];
        return (
            <Cell {...this.props}> 
                <a href={link}>{data[rowIndex][field]}</a>
            </Cell>
       );
    }
}

class TooltipCell extends React.Component {
    render() {
        const {data, rowIndex, columnKey} = this.props;
        //const value = data[rowIndex][columnKey];
        const value = columnKey;
        const toolTip = 'Region: ' + data['region'];
        return (
            <Cell
                {...this.props}
                onMouseEnter={() => { ReactTooltip.show(this.refs.valueDiv); }}
                onMouseLeave={() => { ReactTooltip.hide(this.refs.valueDiv); }}>
                <div ref='valueDiv' data-tip={toolTip}>
                {value}
                </div>
            </Cell>
        );
    }
}


class MyTable extends React.Component {

    constructor(props) {
        super(props);

        this.state = {
            rows: [],
            columns: [],
            filteredEnsembles: [],
            filteredDatasets: [],
            sortBy: 'id',
            sortDir: 'ASC'};
    }

    componentWillMount() {
        this.jsonList();
    }

    jsonList() {
        fetch('/content/ensemble_list').then(
            function(response){
                return response.json();
            }
        ).then(data => {
                var d = data;

                this.setState({
                    rows: d.data,
                    columns: d.columns,
                    filteredEnsembles: d.data,
                    filteredDatasets: d.columns,
                })
            })
    }

    render() {
        var sortDirArrow = '';
        if (this.state.sortDir !== null){
            sortDirArrow = this.state.sortDir === 'DESC' ? ' ↓' : ' ↑';
        }
        var dataset_columns = [];
        for (var index = 0; index < this.state.filteredDatasets.length; index++) {
            var column_tag = this.state.filteredDatasets[index]['dataset'];
            var column_name = column_tag;
            dataset_columns.push(
                <Column 
                    key={column_tag}
                    columnKey={column_tag} 
                    header={<TooltipCell data={this.state.filteredDatasets[index]} columnKey={column_tag}></TooltipCell>}
                    cell={<MyTextCell data={this.state.filteredEnsembles} field={column_tag} />}
                    //cell={<TooltipCell data={this.state.filteredEnsembles} columnKey={column_tag} />}
                    width={120} 
                    allowCellsRecycling={true}
                />
            );
        }
        return (
            <div>
                <br />
                <input
                    onChange={this._onFilterChangeEnsemblesByName.bind(this,'ensemble')}
                    placeholder="Filter Ensembles by Name"
                    style={{width:'200px'}}
                />
                <input 
                    onChange={this._onFilterChangeDatasetsByName.bind(this, 'dataset')}
                    placeholder="Filter Datasets by Name"
                    style={{width:'200px'}}
                />
                <input 
                    onChange={this._onFilterChangeDatasetsByABA.bind(this, 'dataset')}
                    placeholder="Filter Datasets by ABA Region"
                    style={{width:'220px'}}
                />
                <br />
                <Table
                    key={'table'}
                    height={50+((this.state.filteredEnsembles.length+1) * 30)}
                    width={(this.state.filteredDatasets.length * 120) + 250}
                    rowsCount={this.state.filteredEnsembles.length}
                    rowHeight={30}
                    headerHeight={60}>
                    <Column
                        key={'ens_id'}
                        columnKey="id"
                        header={<Cell><a onClick={this._sortRowsBy.bind(this, 'id')}>{'#' + (this.state.sortBy === 'id' ? sortDirArrow: '')}</a></Cell>}
                        cell={<MyTextCell data={this.state.filteredEnsembles} field="id" />}
                        width={50}
                        fixed={true}
                    />
                    <Column 
                        key={'ens_name'}
                        columnKey="ensemble"
                        header={<Cell>{'Ensemble'}</Cell>} 
                        cell={<MyLinkCell data={this.state.filteredEnsembles} field="ensemble" />} 
                        width={200} 
                        fixed={true}
                        //headerRenderer={this._renderHeader.bind(this)}
                    />
                    {dataset_columns}
                </Table>
                <ReactTooltip />
            </div>
        );
    }

    _onFilterChangeEnsemblesByName(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredEnsembles: this.state.rows,
            });
        }else{
            var filterBy = event.target.value.toString().toLowerCase().split(/[ ,]+/);
            var filteredList = [];
            for (var index = 0; index < this.state.rows.length; index++) {
                var v = this.state.rows[index][cellDataKey];
                for (var jndex = 0; jndex < filterBy.length; jndex++) {
                    var compareString = filterBy[jndex];
                    if (v.toString().toLowerCase().indexOf(compareString) !== -1 && compareString !== "") {
                        if (!filteredList.includes(this.state.rows[index])){ 
                            filteredList.push(this.state.rows[index]);
                        }
                    }
                }
            }
            this.setState({
                filteredEnsembles: filteredList,
            });
        }
    }

    _onFilterChangeDatasetsByName(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredEnsembles: this.state.rows,
                filteredDatasets: this.state.columns,
            });
        }else{
            var filterBy = event.target.value.toString().toLowerCase().split(/[ ,]+/);
            var filteredList_ds = [];
            var filteredList_ens = [];
            for (var index_filter = 0; index_filter < filterBy.length; index_filter++) {
                var compareString = filterBy[index_filter];
                for (var index_ds = 0; index_ds < this.state.columns.length; index_ds++) {
                    var v = this.state.columns[index_ds]['dataset'];
                    if (v.toString().toLowerCase().indexOf(compareString) !== -1 && compareString !== "") {
                        if (!filteredList_ds.includes(this.state.columns[index_ds])){ 
                            filteredList_ds.push(this.state.columns[index_ds]);
                        }
                    }
                }
                for (var index_ens = 0; index_ens < this.state.rows.length; index_ens++) {
                    var datasetsInEns = this.state.rows[index_ens]['datasets'].toString().toLowerCase().split(/[\n]+/);

                    for (var jndex = 0; jndex < datasetsInEns.length; jndex++) {
                        if (datasetsInEns[jndex].indexOf(compareString) !== -1 && compareString !== "") {
                            if (!filteredList_ens.includes(this.state.rows[index_ens])) {
                                filteredList_ens.push(this.state.rows[index_ens]);
                            }
                        }
                    }
                }
            }
            this.setState({
                filteredEnsembles: filteredList_ens,
                filteredDatasets: filteredList_ds,
            });
        }

        this.render.bind(this);
    }

    _onFilterChangeDatasetsByABA(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredEnsembles: this.state.rows,
                filteredDatasets: this.state.columns,
            });
        }else{
            var filterBy = event.target.value.toString().toLowerCase().split(/[ ,]+/);
            var size = this.state.columns.length;
            var filteredList_ds = [];

            //Filter through columns (datasets) first
            for (var index_ds = 0; index_ds < size; index_ds++) {
                var v = this.state.columns[index_ds]['region'];
                for (var jndex = 0; jndex < size; jndex++) {
                    var compareString = filterBy[jndex]
                    if (v.toString().toLowerCase().indexOf(compareString) !== -1 && compareString !== "") {
                        if (!filteredList_ds.includes(this.state.columns[index_ds])){ 
                            filteredList_ds.push(this.state.columns[index_ds]);
                        }
                    }
                }
            }
            //Use filtered columns (datasets) list to get list of ensembles (rows) that contain datasets
            var filteredList_ens = [];
            for (var index_ens = 0; index_ens < this.state.rows.length; index_ens++) {
                var datasetsInEns = this.state.rows[index_ens]['datasets'].toString().toLowerCase().split(/[\n]+/);
                for (var index_f = 0; index_f < filteredList_ds.length; index_f++) {
                    var compareString = filteredList_ds[index_f]['dataset'].toString().toLowerCase();
                    if (compareString !== "") {
                        for (var jndex = 0; jndex < datasetsInEns.length; jndex++) {
                            if (datasetsInEns[jndex].indexOf(compareString) !== -1) {
                                if (!filteredList_ens.includes(this.state.rows[index_ens])) {
                                    filteredList_ens.push(this.state.rows[index_ens]);
                                }
                            }
                        }
                    }
                }
            }
            this.setState({
                filteredEnsembles: filteredList_ens,
                filteredDatasets: filteredList_ds,
            });
        }
        this.render.bind(this);
    }

    _sortRowsBy(cellDataKey) {
        var sortDir = this.state.sortDir;
        var sortBy = cellDataKey;
        if (sortBy === this.state.sortBy) {
            sortDir = this.state.sortDir === 'ASC' ? 'DESC' : 'ASC';
        } else {
            sortDir = 'DESC';
        }
        var rows = this.state.filteredEnsembles.slice();
        rows.sort((a, b) => {
            var sortVal = 0;
            if (a[sortBy] > b[sortBy]) {
                sortVal = 1;
            }
            if (a[sortBy] < b[sortBy]) {
                sortVal = -1;
            }

            if (sortDir === 'DESC') {
                sortVal = sortVal * -1;
            }
            return sortVal;
        });

        this.setState({sortBy, sortDir, filteredEnsembles : rows});
    }
}

module.exports = MyTable;
