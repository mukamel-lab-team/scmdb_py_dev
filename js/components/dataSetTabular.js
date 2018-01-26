import React from 'react';
import {Table, Column, Cell} from 'fixed-data-table';
import "fixed-data-table/dist/fixed-data-table.css";

class MyTable extends React.Component {

    constructor() {
        super();

        this.state = {
            rows: [],
            filteredDataList: [],
            sortBy: 'Sample',
            sortDir: 'ASC',
            wheight: '0',
            wwidth: '0'
        };
        this.updateWindowDimensions = this.updateWindowDimensions.bind(this);
    }

    componentWillMount() {
        this.jsonList();
    }
    componentDidMount() {
        this.updateWindowDimensions();
        window.addEventListener('resize', this.updateWindowDimensions);
    }
    componentWillUnmount() {
        window.removeEventListener('resize', this.updateWindowDimensions);
    }
    updateWindowDimensions() {
        this.setState({ wwidth: window.innerWidth, wheight: window.innerHeight });
    }

    jsonList() {
        fetch('/content/metadata/human_MB_EB').then(
            function(response){
                return response.json();
            }
        ).then(data => {
            var d = data;
            this.setState({
                rows: d.data,
                filteredDataList: d.data
            })
        })
    }

    render() {
        var sortDirArrow = '';
        if (this.state.sortDir !== null){
            sortDirArrow = this.state.sortDir === 'DESC' ? ' ↓' : ' ↑';
        }
        return <Table
            height={this.state.wheight - 150}
            width={this.state.wwidth - 100}
            rowsCount={this.state.filteredDataList.length}
            rowHeight={30}
            headerHeight={90}
            rowGetter={function(rowIndex) {return this.state.filteredDataList[rowIndex]; }.bind(this)}>
            <Column dataKey="Sample" width={250} label={'Sample'+ (this.state.sortBy === 'Sample' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Library pool" width={75} label={'Library Pool' + (this.state.sortBy === 'Library pool' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Layer" width={75} label={'Layer' + (this.state.sortBy === 'Layer' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Total reads" width={100} label={'Total Reads' + (this.state.sortBy === 'Total reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Mapped reads" width={100} label={'Mapped reads' + (this.state.sortBy === 'Mapped reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Mapping rate" width={100} label={'Mapping rate' + (this.state.sortBy === 'Mapping rate' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Nonclonal reads" width={100} label={'Nonclonal reads' + (this.state.sortBy === 'Nonclonal reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="% Nonclonal rates" width={150} label={'% Nonclonal rates' + (this.state.sortBy === '% Nonclonal rates' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="mCCC/CCC" width={100} label={'mCCC/CCC' + (this.state.sortBy === 'mCCC/CCC' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="mCG/CG" width={100} label={'mCG/CG' + (this.state.sortBy === 'mCG/CG' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="mCH/CH" width={100} label={'mCH/CH' + (this.state.sortBy === 'mCH/CH' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="Estimated mCG/CG" width={100} label={'Estimated mCG/CG' + (this.state.sortBy === 'Estimated mCG/CG' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="% Genome covered" width={100} label={'% Genome covered' + (this.state.sortBy === '% Genome covered' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="allc file location (Mukamel lab)" width={200} label={'allc file location (Mukamel lab)' + (this.state.sortBy === 'allc file location (Mukamel lab)' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
        </Table>;
    }

    _onFilterChange(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredDataList: this.rows,
            });
        }
        var filterBy = event.target.value.toString().toLowerCase();
        var size = this.rows.length;
        var filteredList = [];
        for (var index = 0; index < size; index++) {
            var v = this.rows[index][cellDataKey];
            if (v.toString().toLowerCase().indexOf(filterBy) !== -1) {
                filteredList.push(this.rows[index]);
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
