import React from 'react';
import {Table, Column, Cell} from 'fixed-data-table';
import "fixed-data-table/dist/fixed-data-table.css";

class MyTable extends React.Component {

    constructor() {
        super();

        this.state = {
            rows: [],
            columns: [],
            columnSize: 0,
            filteredDataList: [],
            sortBy: 'ensemble',
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
                var columnSize = d.columns.length;

                this.setState({
                    rows: d.data,
                    columns: d.columns,
                    filteredDataList: d.data,
                    columnSize
                })
            })
    }

    render() {
        var sortDirArrow = '';
        if (this.state.sortDir !== null){
            sortDirArrow = this.state.sortDir === 'DESC' ? ' ↓' : ' ↑';
        }
        var dataset_columns = []
        for (var index = 0; index < this.state.columnSize; index++) {
            var column_tag = this.state.columns[index]
            var column_name = column_tag
            dataset_columns.push(<Column key={column_tag} dataKey={column_tag} width={100} label={column_name} headerRenderer={this._renderHeader.bind(this)}/>)
        }
        return <Table
            height={52+((this.state.filteredDataList.length+1) * 30)}
            width={(this.state.columnSize * 100) + 200}
            rowsCount={this.state.filteredDataList.length}
            rowHeight={30}
            headerHeight={80}
            rowGetter={function(rowIndex) {return this.state.filteredDataList[rowIndex]; }.bind(this)}>
            <Column dataKey="ensemble" width={200} label={'Ensemble'+ (this.state.sortBy === 'ensemble' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            {dataset_columns}
        </Table>;
    }
    _onFilterChange(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredDataList: this.state.rows,
            });
        }
        var filterBy = event.target.value.toString().toLowerCase();
        var size = this.state.rows.length;
        var filteredList = [];
        for (var index = 0; index < size; index++) {
            var v = this.state.rows[index][cellDataKey];
            if (v.toString().toLowerCase().indexOf(filterBy) !== -1) {
                filteredList.push(this.state.rows[index]);
            }
        }
        this.setState({
            filteredDataList: filteredList,
        });
    }
    _sortRowsBy(cellDataKey) {
        var sortDir = this.state.sortDir;
        var sortBy = cellDataKey;
        if (sortBy === this.state.sortBy) {
            sortDir = this.state.sortDir === 'ASC' ? 'DESC' : 'ASC';
        } else {
            sortDir = 'DESC';
        }
        var rows = this.state.filteredDataList.slice();
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

        this.setState({sortBy, sortDir, filteredDataList : rows});
    }

    _renderHeader(label, cellDataKey) {
        return <div>
            <a onClick={this._sortRowsBy.bind(this, cellDataKey)}>{label}</a>
            <div>
                <input type="text" style={{width:90+'%'}} onChange={this._onFilterChange.bind(this, cellDataKey)}/>
            </div>
        </div>;
    }
}

module.exports = MyTable;
